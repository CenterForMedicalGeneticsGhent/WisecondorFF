import copy
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Generator, Optional, Tuple

import numpy as np
import pysam


def read_pair_gen(
    bam_pointer: pysam.AlignmentFile,
    stats_map: Dict[str, Any],
    region_str: Optional[str] = None,
) -> Generator:
    """
    Return read pairs as they are encountered,
    includes the position of the first previously aligned reads
    """

    read_dict = defaultdict(lambda: [None, None])
    lr_pos = -1

    for read in bam_pointer.fetch(region=region_str):
        stats_map["reads_seen"] += 1

        if not read.is_proper_pair:
            stats_map["reads_pairf"] += 1
            continue

        if read.is_secondary or read.is_supplementary:
            stats_map["reads_secondary"] += 1
            continue

        query_name = read.query_name
        if query_name not in read_dict:
            if read.is_read1:
                read_dict.get(query_name)[0] = (read, lr_pos)
            else:
                read_dict.get(query_name)[1] = (read, lr_pos)
        else:
            if read.is_read1:
                yield (read, lr_pos), read_dict[query_name][1]
            else:
                yield read_dict[query_name][0], (read, lr_pos)
            del read_dict[query_name]
        lr_pos = read.pos


def convert_reads(
    reads_path: Path,
    binsize: int = 5000,
    min_mapq: int = 1,
    reference_fasta: Optional[Path] = None,
) -> Tuple[Dict[str, Dict[Any, Any]], Dict[str, int]]:
    """
    Convert reads to a dictionary of samples, each containing a dictionary of chromosomes,
    """
    stats_map = defaultdict(int)
    chr_re = re.compile("chr", re.IGNORECASE)

    rc_chr_map = {}
    for chrom in map(str, range(1, 23)):
        rc_chr_map[chrom] = None
    fs_chr_map = copy.deepcopy(rc_chr_map)

    def get_midpoint(r1, r2):
        return ((r1.reference_start + r2.reference_end) // 2), r1.template_length

    if reads_path.suffix() == ".bam":
        reads_file = pysam.AlignmentFile(str(reads_path), "rb")
    elif reads_path.suffix() == ".cram":
        if reference_fasta is not None:
            reads_file = pysam.AlignmentFile(
                str(reads_path), "rc", reference_filename=str(reference_fasta)
            )
        else:
            logging.error(
                "Cram support requires a reference file, please use the --reference-fasta argument"
            )
            sys.exit(1)
    else:
        logging.error(
            "Unsupported input file type. Make sure your input filename has a correct extension ( bam or cram)"
        )
        sys.exit(1)

    with reads_file:
        try:
            reads_file.check_index()
        except ValueError:
            logging.info("No index found! Trying to generate...")
            pysam.index(str(reads_path))
            sys.exit(2)
        except AttributeError:
            logging.info("File is SAM formatted and thus has no index.")
            sys.exit(3)

        for index, chrom in enumerate(reads_file.references):
            chrom_name = chr_re.sub("", chrom)
            if chrom_name not in rc_chr_map:
                continue

            chr_size = reads_file.lengths[index]
            n_bins = (int(chr_size / float(binsize))) + 1

            logging.info(f"Processing: {chrom}, filling: {n_bins} bins")

            fs_dict = [defaultdict(int) for __ in range(n_bins)]
            rc_counts = np.zeros(
                int(reads_file.lengths[index] / float(binsize) + 1), dtype=np.int32
            )

            for ((read1, read1_prevpos), (read2, read2_prevpos)) in read_pair_gen(
                reads_file, stats_map, chrom
            ):
                if read1.pos == read1_prevpos or read2.pos == read2_prevpos:
                    stats_map["reads_rmdup"] += 2
                    continue

                if read1.mapping_quality < min_mapq or read2.mapping_quality < min_mapq:
                    stats_map["reads_mapq"] += 2
                    continue

                stats_map["reads_kept"] += 2

                rc_counts[int(read1.pos / binsize)] += 1
                rc_counts[int(read2.pos / binsize)] += 1

                if read1.template_length > 0:
                    mid, insert_size = get_midpoint(read1, read2)
                else:
                    mid, insert_size = get_midpoint(read2, read1)
                fs_dict[int(mid // binsize)][insert_size] += 1

            rc_chr_map[chrom_name] = rc_counts
            fs_chr_map[chrom_name] = fs_dict

        qual_info = {
            "mapped": reads_file.mapped,
            "unmapped": reads_file.unmapped,
            "no_coordinate": reads_file.nocoordinate,
            "filter_rmdup": stats_map["reads_rmdup"],
            "filter_mapq": stats_map["reads_mapq"],
            "pre_retro": stats_map["reads_seen"],
            "post_retro": stats_map["reads_kept"],
            "pair_fail": stats_map["reads_pairf"],
        }

    return {"RC": rc_chr_map, "FS": fs_chr_map}, qual_info
