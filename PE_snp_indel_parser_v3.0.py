#!/usr/bin/env python3
## version 3.0 (05/23/2025, Simone Gupta)
## NM tag < 20
## the mistamtch error due to non integers. 
## bugs corrected (05/23/2025, Simone Gupta)
## Insertion report <ref_pos>:<length>:<average base quality>:<average phred error>
## added Process only Primary alignment 
## Deduplicating insertions same reference position multiple times. This can happen in complex CIGAR strings (like 3I6M1D4M2I5M) when insertions are adjacent or overlap in the read, but the reference position is not incremented after an insertion (I doesn't move the reference position)



import pysam
import pandas as pd
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed



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

            avg_qual = round(sum(qual) / len(qual), 2) if qual else None
            avg_errP = round(sum(error_P) / len(error_P), 5) if error_P else None

            indels_mismatches['I'].append((ref_pos, length, avg_qual, avg_errP))
            read_edit_quality['I'].append(avg_qual)
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
        seen_positions = set()

        for pos, length, qual, error_P in indels[type_]:
            if pos in seen_positions:
                continue
            seen_positions.add(pos)

            formatted[type_].append(f"{pos}:{length}:{qual}:{error_P}")

        formatted[type_] = ";".join(formatted[type_])
    return formatted

def format_qcindel(qc_indel):
    return qc_indel['I'][0], qc_indel['D'][0]

def process_reads(bamfile_path, sample, chunk):
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    rows = []
    wt_reads = 0
    read_pair_count = 0

    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        
        tags = dict(read.tags)
        if tags.get('tp') != 'P' or read.is_unmapped or read.query_sequence is None:
            continue
        
        if tags.get('NM', 0) > 20:
            continue

        if read.query_sequence.count("N") > 3 or read.mapping_quality < 30:
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
        read_pair_count += 1

    bamfile.close()
    result_df = pd.DataFrame(rows)

    result_df['Mismatches'] = result_df['Mismatches'].fillna('').astype(str).str.strip()

    excluded_mask = (
        (result_df['INS_QC'] == 'I_NA') & 
        (result_df['DEL_QC'] == 'D_NA') & 
        (result_df['Mismatches'] == '')
    )
    excluded_count = excluded_mask.sum()
    wt_reads += excluded_count

    result_df = result_df[~excluded_mask]
    result_df['WildCount'] = wt_reads
    result_df['AlignedPairs'] = read_pair_count

    return result_df

def chunk_positions(values, chunk_size):
    flat_positions = [pos for sublist in values for pos in sublist]
    for i in range(0, len(flat_positions), chunk_size):
        yield flat_positions[i:i + chunk_size]

def main():
    bamfile_path = sys.argv[1]  # Path to BAM file
    sample = sys.argv[2]

    
    # Split BAM based on depth, now we get chunks containing reads and mate indices
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    dict_reads = defaultdict(list)
  

    # Fetch reads and store their indices
    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        dict_reads[read.query_name].append(idx)
    read_pair_count = len(dict_reads)

    chunks = list(chunk_positions(list(dict_reads.values()), 2500))

    num_processes = 8  # Adjust based on your CPU
    
    wt_reads = 0
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_chunk = {executor.submit(process_reads, bamfile_path, sample, chunk): (chunk) 
                           for chunk in (chunks)}
        
        # Collect the results
        all_dfs = []
        for future in as_completed(future_to_chunk):
            #print(future.result())
            # Group by 'Query_Name' and apply the concatenation function
            result_df = future.result()
            result_df['Mismatches'] = result_df['Mismatches'].fillna('').astype(str).str.strip()

            # Filter the DataFrame
            excluded_mask = (
            (result_df['INS_QC'] == 'I_NA') & 
            (result_df['DEL_QC'] == 'D_NA') & 
            (result_df['Mismatches'] == '')
            )
            excluded_count = excluded_mask.sum()  # Number of rows excluded

            result_df = result_df[~excluded_mask]
            all_dfs.append(result_df)
            wt_reads = wt_reads +  excluded_count
    # Concatenate all DataFrames
    final_df = pd.concat(all_dfs, ignore_index=True)
    final_df['WildCount'] = wt_reads
    final_df['AlignedPairs'] = read_pair_count

    # Display the resulting DataFrame
    # Save the result to a file
    outfile = f"{sample}_snp_indel.output_data.tsv"
    out_dir = os.getcwd()
    outfile = os.path.join(out_dir, outfile)
    final_df.to_csv(outfile, sep='\t', index=False)
    print(f"Output saved to {outfile}")
if __name__ == "__main__":
    main()
