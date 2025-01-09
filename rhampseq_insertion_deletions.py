import pysam
import sys
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from array import array
from collections import defaultdict

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


def extract_indels_and_mismatches(read):
    """
    Extract multiple insertions (I), deletions (D), and mismatches (X) from CIGAR and MD tags.
    Handles reverse strand alignment by adjusting coordinates.
    """
    # Initialize output dictionary
    indels_mismatches = {'I': [], 'D': [], 'X': []}
    
    # Start position and quality scores
    ref_pos = read.reference_start
    base_qualities = read.query_qualities
    read_pos = 0

    # Parse CIGAR string
    for cigar_op, length in read.cigar:
        if cigar_op == 0:  # Match (M) or mismatch
            ref_pos += length
            read_pos += length
        elif cigar_op == 1:  # Insertion (I)
            # Capture insertion
            quality = base_qualities[read_pos:read_pos + length]
            error_P = [10 ** (-q / 10) for q in quality]  # Probability of error
            indels_mismatches['I'].append((ref_pos, length, quality, error_P))
            read_pos += length
        elif cigar_op == 2:  # Deletion (D)
            # Capture deletion
            quality = base_qualities[read_pos - 1] if read_pos > 0 else None
            error_P = 10 ** (-quality / 10) if quality is not None else None
            indels_mismatches['D'].append((ref_pos, length, quality, error_P))
            ref_pos += length  # Move reference position forward

    # Parse MD tag for mismatches (X)
    md = dict(read.tags).get('MD')
    if md:
        ref_offset = 0
        num = ''
        for char in md:
            if char.isdigit():
                num += char  # Accumulate numbers (matches)
            else:
                # Move reference position forward for matches
                if num:
                    ref_offset += int(num)
                    num = ''
                # Handle mismatches (X)
                if char.isalpha():  # Mismatch
                    mismatch_pos = ref_pos + ref_offset
                    if read_pos < len(base_qualities):  # Boundary check
                        quality = base_qualities[read_pos]
                        error_P = 10 ** (-quality / 10)
                    else:
                        quality = "NA"
                        error_P = "NA"
                    indels_mismatches['X'].append((mismatch_pos, 1, quality, error_P))
                    ref_offset += 1
                    read_pos += 1  # Update read position
                elif char == '^':  # Deletion marker
                    num = ''  # Skip deletions handled in CIGAR
    #print(indels_mismatches)
    return indels_mismatches

def hamming(x, y):
    #print(f"Hamming: {x}, {y}")  # Print the sequences being compared
    return sum(a == b for a, b in zip(x, y)) / len(x)   ## similarity ratio (proportion of matching elements) between two sequences x and y

def merge_to_bam(read_sequence, mate_sequence, min_overlap=2, quality_score=30):
    # Handle read orientation. I am not using the quality score for now

    seq1, seq2 = read_sequence, mate_sequence
    max_len = len(seq1)

    # Merge sequences based on overlap
    for i in range(max_len - min_overlap + 1):
        overlap_len = max_len - i
        if overlap_len < min_overlap: break
        identity = hamming(seq1[i:], seq2[:overlap_len])
        if identity >= 0.5:
            merged_seq = seq1[:i] + ''.join(
                s1 if s1 == s2 else s2 for s1, s2 in zip(seq1[i:], seq2[:overlap_len])
            ) + seq2[overlap_len:]
             
            return overlap_len, 1  # Return 1 indicating that the merge was successful

    # Return 0 if no merge was found
    return 0, 0


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

def process_reads(bamfile_path, sample, contig, start, end,chunk):
    
    """
    Process reads in the given range (start, end) and extract relevant information.
    This function is called in parallel for each chunk.
    """
    # Open the BAM file using pysam
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")

    # Initialize DataFrame for results
    columns = [
        "Query_Name", "ReadReferenceId", "Reference_Start", "Reference_End", "CIGAR_String", "Orientation",
        "QC_Fail", "Secondary", "Supplementary", "Duplicate", "Insertions", "Deletions", "Mismatches", "Sequence", "Read_1", "Read_2"
       
    ]
    df = pd.DataFrame(columns=columns)  # Initialize an empty DataFrame

    
    
    # Fetch reads from the specified region
    for idx, read in enumerate(bamfile.fetch(contig, start, end, until_eof=True)):
        for position_index in chunk:
        
            if idx == position_index:
                #print(f"Read {idx}: {read}")
        
                if read.is_unmapped:
                    continue
                if read.mapping_quality >= 30:
                    is_reverse = read.is_reverse
                    orientation = "reverse" if is_reverse else "forward"
                    # Extract indels and mismatches for both read and mate
                    indels_mismatches = extract_indels_and_mismatches(read)
                    # Merge reads if needed (Optional, based on your logic)
                    #overlap_len, merge_status = merge_to_bam(read, mate)
                    #merge_overlap = overlap_len
                
                    # Format indels/mismatches into columns
                    read_indels = format_indels(indels_mismatches)
                    # Write details for the read and its mate in one row
                    infor_row = {
                    "Query_Name": read.query_name,
                    "ReadReferenceId": read.reference_name,
                    "Reference_Start": read.reference_start,
                    "Reference_End": read.reference_end,
                    "CIGAR_String": read.cigarstring,
                    "Orientation": orientation,
                    "QC_Fail": read.is_qcfail,
                    "Secondary": read.is_secondary,
                    "Supplementary": read.is_supplementary,
                    "Duplicate": read.is_duplicate,
                    "Insertions": read_indels['I'],
                    "Deletions": read_indels['D'],
                    "Mismatches": read_indels['X'],
                    "Sequence": read.query_sequence,
                    "Read_1": read.is_read1,
                    "Read_2": read.is_read2
                    }
                    df = pd.concat([df, pd.DataFrame([infor_row])], ignore_index=True)

    bamfile.close()
    return df

def process_merge_for_reads(df):
    # Initialize lists to hold the overlap lengths and merge statuses
    overlap_lens = []
    merge_successes = []

    # Iterate through each row of the dataframe
    for idx, merged_row in df.iterrows():
        # Check if both Read_Sequence and Mate_Sequence are not 'NA'
        read_sequence = merged_row.get("Read_Sequence", 'NA')
        mate_sequence = merged_row.get("Mate_Sequence", 'NA')

        # Only proceed with the merge if both sequences are not 'NA'
        if read_sequence != 'NA' and mate_sequence != 'NA':
            overlap_len, merge_status = merge_to_bam(read_sequence, mate_sequence, min_overlap=2, quality_score=30)
        else:
            overlap_len = 0
            merge_status = 0  # No merge if one or both sequences are 'NA'

        # Append the results to the lists
        overlap_lens.append(overlap_len)
        merge_successes.append(merge_status)

    # Add the results as new columns in the dataframe
    df['Overlap_Len'] = overlap_lens
    df['Merge_Success'] = merge_successes

    return df

# Function to restructure the DataFrame based on the provided criteria
def restructure_paired_end(df):
    # Split the data into Read_1 (forward) and Read_2 (reverse) based on Flag
    read_1_df = df[df['Read_1'] == True]
    read_2_df = df[df['Read_2'] == True]
    
    # Initialize a list to hold the final rows
    final_rows = []
    
    # Loop through each unique Query_Name and merge corresponding rows
    for query_name in df['Query_Name'].unique():
        # Get the Read_1 (first) and Read_2 (second) rows for the given Query_Name
        read_1_row = read_1_df[read_1_df['Query_Name'] == query_name]
        read_2_row = read_2_df[read_2_df['Query_Name'] == query_name]
        
        # Create a dictionary to hold the merged row
        merged_row = {}
        
        # For Read_1 columns (columns without 'Mate_' prefix)
        for col in read_1_row.columns:
            if col not in ['Read_1', 'Read_2']:
                merged_row[f"Read_{col}"] = read_1_row[col].values[0] if not read_1_row.empty else 'NA'
        
        # For Read_2 columns (columns with 'Mate_' prefix)
        for col in read_2_row.columns:
            if col not in ['Read_1', 'Read_2']:
                merged_row[f"Mate_{col}"] = read_2_row[col].values[0] if not read_2_row.empty else 'NA'
        # Add the merged row to the final rows list
        final_rows.append(merged_row)
    
    # Create the final DataFrame
    result_df = pd.DataFrame(final_rows)
    
    return result_df



def chunk_positions(values, chunk_size):
    # Flatten the list of positions (assuming each value is a list of positions)
    flat_positions = [pos for sublist in values for pos in sublist]
    
    # Split the flat list into chunks
    for i in range(0, len(flat_positions), chunk_size):
        #print(f"Starting index for chunk: {i}")
        yield flat_positions[i:i + chunk_size]

def main():
    bam_location = sys.argv[1]  # Path to BAM file
    sample = sys.argv[2]
    contig = sys.argv[3]
    start = int(sys.argv[4])
    end = int(sys.argv[5])
    
    bamfile_path = f"{bam_location}/sorted_{sample}_aligned.bam"
    # Split BAM based on depth, now we get chunks containing reads and mate indices
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    dict_reads = defaultdict(list)
    
    # Fetch reads and store their indices
    for idx, read in enumerate(bamfile.fetch(contig, start, end, until_eof=True)):
        dict_reads[read.query_name].append(idx)


    #Now we have the reads, mate pair and their indices, we can chunk them
    chunks = chunk_positions(list(dict_reads.values()), 2500)

    # Print or process each chunk
    #for chunk in chunks:
        #print(f"Processing chunk: {chunk}")

    num_processes = 8  # Adjust based on your CPU
    

    # Process_reads in Chunks
    #print(f"Processing reads from {chunk_starts} to {chunk_ends}")
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_chunk = {executor.submit(process_reads, bamfile_path, sample, contig, start, end, chunk): (chunk) 
                           for chunk in (chunks)}


        
        # Collect the results
        all_dfs = []
        for future in as_completed(future_to_chunk):
            #print(future.result())
            # Group by 'Query_Name' and apply the concatenation function
            result_df = future.result()
            result_df = restructure_paired_end(result_df)
            # add the merge status
            result_df = process_merge_for_reads(result_df)
            all_dfs.append(result_df)

    # Concatenate all DataFrames
    final_df = pd.concat(all_dfs, ignore_index=True)
    

    # Display the resulting DataFrame
    # Save the result to a file
    outfile = f"{sample}_{contig}_{start}_{end}.output_data.tsv"
    out_dir = "/lrlhps/genomics/prod/lgm/dna_editing/TTR_rhAmp_LOD/BN24-8124/scripts/"
    outfile = out_dir + outfile
    final_df.to_csv(outfile, sep='\t', index=False)
    print(f"Output saved to {outfile}")
if __name__ == "__main__":
    main()
