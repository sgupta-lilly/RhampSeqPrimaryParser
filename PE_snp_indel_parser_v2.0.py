#!/usr/bin/env python3
## version 2.0 (05/23/2025, Simone Gupta)
## Insertion report <ref_pos>:<length>:<average base quality>:<average phred error>
## added Process only Primary alignment 
## Deal with deduplicating insertions same reference position multiple times. This can happen in complex CIGAR strings (like 3I6M1D4M2I5M) when insertions are adjacent or overlap in the read, but the reference position is not incremented after an insertion (I doesn't move the reference position)

import pysam
import sys
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import re

def count_mismatches(read):
    md = dict(read.tags).get('MD')
    return sum(1 for char in md if char.isalpha() and char not in '^') if md else 0

def reference_quality(read):
    query_dict_quality, query_dict_reference = [], []
    for query_idx, ref_idx in read.get_aligned_pairs():
        if ref_idx is None:
            continue
        if query_idx is None:
            query_dict_quality.append(0)
            query_dict_reference.append(ref_idx)
        else:
            q = read.query_qualities[query_idx] if read.query_qualities else None
            query_dict_quality.append(q)
            query_dict_reference.append(ref_idx)
    return query_dict_quality, query_dict_reference

def classify_read(read_edit_quality, quality_threshold=29):
    classify_IDX = {'I': [], 'D': []}
    for key in ['I', 'D']:
        if not read_edit_quality[key]:
            classify_IDX[key].append(f"{key}_NA")
        else:
            qual = read_edit_quality[key][0]
            if qual is None:
                classify_IDX[key].append(f"{key}_NA")
            elif qual > quality_threshold:
                classify_IDX[key].append(f"{key}_High")
            else:
                classify_IDX[key].append(f"{key}_Low")
    return classify_IDX    

def extract_indels_and_mismatches(read):
    indels_mismatches = {'I': [], 'D': [], 'X': []}
    read_edit_quality = {'I': [], 'D': [], 'X': []}
    ref_pos, read_pos = read.reference_start, 0
    base_qualities = read.query_qualities if read.query_qualities is not None else []

    for cigar_op, length in read.cigar:
        if cigar_op == 4:  # Soft clip
            read_pos += length
        elif cigar_op == 0:  # Match
            ref_pos += length
            read_pos += length
        elif cigar_op == 1:  # Insertion
            qual = base_qualities[read_pos:read_pos + length] if base_qualities else [None] * length
            error_P = [10 ** (-q / 10) if q is not None else None for q in qual]
            indels_mismatches['I'].append((ref_pos, length, qual, error_P))
            read_edit_quality['I'].append(min(qual) if qual else None)
            read_pos += length
        elif cigar_op == 2:  # Deletion
            qual = base_qualities[read_pos - 1] if read_pos > 0 and base_qualities else None
            error_P = 10 ** (-qual / 10) if qual is not None else None
            indels_mismatches['D'].append((ref_pos, length, qual, error_P))
            read_edit_quality['D'].append(qual)
            ref_pos += length

    # MD tag for mismatches
    try:
        md_tag = read.get_tag('MD')
    except KeyError:
        md_tag = ""
    ref_pos, read_pos = read.reference_start, 0
    md_tokens = re.findall(r'(\d+|\^[A-Z]+|[A-Z])', md_tag)
    for token in md_tokens:
        if token.isdigit():
            n = int(token)
            ref_pos += n
            read_pos += n
        elif token.startswith('^'):
            ref_pos += len(token) - 1
        else:
            qual = base_qualities[read_pos] if read_pos < len(base_qualities) else None
            error_P = 10 ** (-qual / 10) if qual is not None else None
            indels_mismatches['X'].append((ref_pos, 1, qual, error_P))
            read_edit_quality['X'].append(qual)
            ref_pos += 1
            read_pos += 1

    QC_INDEL = classify_read(read_edit_quality)
    return indels_mismatches, QC_INDEL

def format_indels(indels):
    formatted = {}
    for type_ in ['I', 'D', 'X']:
        formatted[type_] = []
        seen = set()
        for pos, length, quality, error_P in indels[type_]:
            if type_ == 'I':
                for i in range(length):
                    q = quality[i] if i < len(quality) else None
                    ep = error_P[i] if i < len(error_P) else None
                    entry = f"{pos + i}:1:{q}:{ep}"
                    if entry not in seen:
                        seen.add(entry)
                        formatted[type_].append(entry)
            else:
                entry = f"{pos}:{length}:{quality}:{error_P}"
                if entry not in seen:
                    seen.add(entry)
                    formatted[type_].append(entry)
        formatted[type_] = ";".join(formatted[type_])
    return formatted

def format_qcindel(qc_indel):
    return qc_indel['I'][0], qc_indel['D'][0]

def process_reads(bamfile_path, sample, chunk):
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    rows = []

    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        tags = dict(read.tags)
        if tags.get('tp') != 'P' or read.is_unmapped or read.query_sequence is None:
            continue
        if read.query_sequence.count("N") > 3 or read.mapping_quality < 30:
            continue
        if idx not in chunk:
            continue

        orientation = "reverse" if read.is_reverse else "forward"
        indels_mismatches, qc_indel = extract_indels_and_mismatches(read)
        INS_QC, DEL_QC = format_qcindel(qc_indel)
        read_indels = format_indels(indels_mismatches)
        read_quality, read_reference = reference_quality(read)

        rows.append({
            "Query_Name": read.query_name,
            "ReadReferenceId": read.reference_name,
            "Reference_Start": read.reference_start,
            "Reference_End": read.reference_end,
            "CIGAR_String": read.cigarstring,
            "Orientation": orientation,
            "Insertions": read_indels['I'],
            "Deletions": read_indels['D'],
            "Mismatches": read_indels['X'],
            "Sequence": read.query_sequence,
            "Read_1": read.is_read1,
            "Read_2": read.is_read2,
            "INS_QC": INS_QC,
            "DEL_QC": DEL_QC,
            "ReadQualities": read_quality,
            "ReferencePosition": read_reference
        })
    bamfile.close()
    return pd.DataFrame(rows)

def chunk_positions(values, chunk_size):
    flat_positions = [pos for sublist in values for pos in sublist]
    for i in range(0, len(flat_positions), chunk_size):
        yield flat_positions[i:i + chunk_size]

def main():
    bamfile_path = sys.argv[1]
    sample = sys.argv[2]

    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    dict_reads = defaultdict(list)
    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        dict_reads[read.query_name].append(idx)
    bamfile.close()

    chunks = list(chunk_positions(list(dict_reads.values()), 2500))
    num_processes = 8

    final_df = pd.DataFrame()
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_reads, bamfile_path, sample, chunk): chunk for chunk in chunks}
        for future in as_completed(futures):
            df_chunk = future.result()
            final_df = pd.concat([final_df, df_chunk], ignore_index=True)

    output_file = f"{sample}_read_summary.csv"
    final_df.to_csv(output_file, index=False)
    print(f"Finished. Output saved to: {output_file}")

if __name__ == "__main__":
    main()
