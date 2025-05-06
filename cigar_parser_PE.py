import pysam
import sys
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from array import array
from collections import defaultdict
import re
def count_mismatches(read):
    """
    Count the number of mismatches from the MD tag.
    Mismatches are counted from alphabetic characters in the MD tag (A, C, G, T).
    Deletions in the MD tag denoted by '^'are not counted as mismatches.
    """
    md = dict(read.tags).get('MD')  # Get the MD tag
    mismatch_count = 0  # Initialize mismatch count

    # Count mismatches from MD tag
    if md:
        # Mismatches are represented by alphabetic characters in the MD tag (A, C, G, T)
        mismatch_count = sum(1 for char in md if char.isalpha() and char not in '^')

    return mismatch_count


# Fetch read quality and reference sequence
def reference_quality(read):
    query_dict_quality = []
    query_dict_reference = []

    for each in read.get_aligned_pairs():
        if each[1] is None:
            continue
        if each[0] is None:
            query_dict_quality.append(0)
            query_dict_reference.append(each[1])
        elif each[0] is not None and each[1] is not None:
            # Check if query_qualities is None
            if read.query_qualities is not None:
                query_dict_quality.append(read.query_qualities[each[0]])
            else:
                query_dict_quality.append(None)  # Append None if no quality scores are available
            query_dict_reference.append(each[1])

    return query_dict_quality, query_dict_reference

def classify_read(read_edit_quality, quality_threshold=29):
    classify_IDX = {'I': [], 'D': []}
    
    for key in read_edit_quality:
        if key == 'I':
            if read_edit_quality[key] == []:
                classify_IDX[key].append("I_NA")
            if read_edit_quality[key] != []:
                if (read_edit_quality[key][0]) > quality_threshold:
                   classify_IDX[key].append("I_High")
                elif (read_edit_quality[key][0]) <= quality_threshold:   
                    classify_IDX[key].append("I_Low")
        if key == 'D':
            if read_edit_quality[key] == []:
                classify_IDX[key].append("D_NA")
            if read_edit_quality[key] != []:
                if (read_edit_quality[key][0]) > quality_threshold:
                   classify_IDX[key].append("D_High")
                elif (read_edit_quality[key][0]) <= quality_threshold:   
                    classify_IDX[key].append("D_Low")
    return classify_IDX    
import re

def extract_indels_and_mismatches(read):
    """
    Extracts insertions (I), deletions (D), and mismatches (X) from CIGAR and MD tags.
    Handles read positions and reference positions correctly.
    Returns:
      - indels_mismatches: dictionary of { 'I': [...], 'D': [...], 'X': [...] }
      - read_edit_quality: quality summary used for QC
    """
    indels_mismatches = {'I': [], 'D': [], 'X': []}
    read_edit_quality = {'I': [], 'D': [], 'X': []}

    ref_pos = read.reference_start
    read_pos = 0
    base_qualities = read.query_qualities if read.query_qualities is not None else []

    # First parse CIGAR for insertions (I) and deletions (D)
    for cigar_op, length in read.cigar:
        if cigar_op == 4:  # Soft clip (S)
            read_pos += length
        elif cigar_op == 0:  # Match (M)
            ref_pos += length
            read_pos += length
        elif cigar_op == 1:  # Insertion (I)
            qual = base_qualities[read_pos:read_pos + length] if base_qualities else [None] * length
            error_P = [10 ** (-q / 10) if q is not None else None for q in qual]
            indels_mismatches['I'].append((ref_pos, length, qual, error_P))
            read_edit_quality['I'].append(min(qual) if qual else None)
            read_pos += length
        elif cigar_op == 2:  # Deletion (D)
            qual = base_qualities[read_pos - 1] if read_pos > 0 and base_qualities else None
            error_P = 10 ** (-qual / 10) if qual is not None else None
            indels_mismatches['D'].append((ref_pos, length, qual, error_P))
            read_edit_quality['D'].append(qual)
            ref_pos += length

    # Then parse MD tag for mismatches (X)
    try:
        md_tag = read.get_tag('MD')
    except KeyError:
        md_tag = ""

    ref_pos = read.reference_start
    read_pos = 0
    md_tokens = re.findall(r'(\d+|\^[A-Z]+|[A-Z])', md_tag)

    for token in md_tokens:
        if token.isdigit():
            matches = int(token)
            ref_pos += matches
            read_pos += matches
        elif token.startswith('^'):
            deleted_bases = token[1:]
            ref_pos += len(deleted_bases)
        else:
            # Mismatch
            if read_pos < len(base_qualities):
                qual = base_qualities[read_pos]
                error_P = 10 ** (-qual / 10)
            else:
                qual = None
                error_P = None
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
        
        for pos, length, quality, error_P in indels[type_]:
            # Convert quality and error_P to lists if they are not already
            quality = list(quality) if isinstance(quality, (array, bytes)) else quality
            error_P = list(error_P) if isinstance(error_P, (array, bytes)) else error_P
            
            # Debugging: Check if quality and error_P are lists and print them
            #print(f"Processing {type_} at position {pos} with length {length}...")
            #print(f"Quality: {quality}, Error_P: {error_P}")
            
            if type_ == 'I' and isinstance(quality, list) and isinstance(error_P, list):
                # Handle multiple base insertions where quality and error_P are lists
                for i in range(length):
                    formatted[type_].append(f"{pos + i}:{1}:{quality[i]}:{error_P[i]}")
            else:
                # Handle cases where quality and error_P are single values or arrays of length 1
                formatted[type_].append(f"{pos}:{length}:{quality}:{error_P}")
        
        # Join the formatted values with semicolons
        formatted[type_] = ";".join(formatted[type_])
    
    return formatted
def format_qcindel(qc_indel):
    INS_QC = qc_indel['I'][0]
    DEL_QC = qc_indel['D'][0]
    return INS_QC, DEL_QC

def process_reads(bamfile_path, sample, chunk):
    
    """
    Process reads in the given range (start, end) and extract relevant information.
    This function is called in parallel for each chunk.
    """
    # Open the BAM file using pysam
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")

    # Initialize DataFrame for results
    columns = [
        "Query_Name", "ReadReferenceId", "Reference_Start", "Reference_End", "CIGAR_String", "Orientation",
         "Insertions", "Deletions", "Mismatches", "Sequence", "INS_QC", "DEL_QC", "ReadQualities", "ReferencePosition"
       
    ]
    df = pd.DataFrame(columns=columns)  # Initialize an empty DataFrame

    
    
    # Fetch reads from the specified region
    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        for position_index in chunk:
        
            if idx == position_index:
                #print(f"Read {idx}: {read}")
        
                if read.is_unmapped:
                    continue
                if read.query_sequence.count("N") > 3:
                    continue
                if read.mapping_quality >= 30:
                    is_reverse = read.is_reverse
                    orientation = "reverse" if is_reverse else "forward"
                    # Extract indels and mismatches for both read and mate
                    indels_mismatches,qc_indel = extract_indels_and_mismatches(read)
                    INS_QC, DEL_QC = format_qcindel(qc_indel)
                
                    # Format indels/mismatches into columns
                    read_indels = format_indels(indels_mismatches)
                    read_quality, read_reference = reference_quality(read)
                    # Write details for the read and its mate in one row
                    infor_row = {
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
                    "INS_QC": INS_QC,
                    "DEL_QC": DEL_QC,
                    "ReadQualities": read_quality,
                    "ReferencePosition": read_reference
                    }
                    df = pd.concat([df, pd.DataFrame([infor_row])], ignore_index=True)

    bamfile.close()
    return df

def chunk_positions(values, chunk_size):
    # Flatten the list of positions (assuming each value is a list of positions)
    flat_positions = [pos for sublist in values for pos in sublist]
    
    # Split the flat list into chunks
    for i in range(0, len(flat_positions), chunk_size):
        #print(f"Starting index for chunk: {i}")
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

    #Create Chunks of ~5000 reads
    #Now we have the reads, mate pair and their indices, we can chunk them
    chunks = list(chunk_positions(list(dict_reads.values()), 2500))

    # Print or process each chunk
    #for chunk in chunks:
        #print(f"Processing chunk: {chunk}")

    num_processes = 8  # Adjust based on your CPU
    
    wt_reads = 0
    # Process_reads in Chunks
    #print(f"Processing reads from {chunk_starts} to {chunk_ends}")
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

            # Count number of rows that meet the exclusion filter


            # Add excluded_count as a column to all rows
            #result_df['Excluded_Count'] = excluded_count
   
            #result_df = restructure_paired_end(result_df)
            # add the merge status
            #result_df = process_merge_for_reads(result_df)
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
