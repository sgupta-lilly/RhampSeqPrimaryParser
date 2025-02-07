import sys
import os
import pysam
import pandas as pd
import re  # Needed for parsing CIGAR string

def write_chunk_to_file(results, output_file, write_header=False):
    """Writes a chunk of data to a TSV file, ensuring headers are written only once."""
    df = pd.DataFrame(results, columns=[
        "Read Name", "Read Chr", "Read Start", "Mate Chr", "Mate Start", 
        "Read CIGAR", "Mate CIGAR", "Read MAPQ", "Mate MAPQ", 
        "Read AS", "Mate AS", "Event Type", "Read SA", "Mate SA"])
    
    df.to_csv(output_file, sep="\t", index=False, mode="a", header=write_header)
    results.clear()  # Free memory

def extract_high_confidence_chimeric_reads(bamfile_path, output_file, chunk_size=10000, distance_threshold=1000, deletion_threshold=50, mapq_threshold=30, alignment_score_threshold=100):
    """Extracts chimeric reads using chunked processing to optimize memory usage and detect long deletions."""
    
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    
    read_pairs = {}  # Dictionary to store read-mate pairs temporarily
    results = []  # List to store results before writing in chunks

    # Remove existing file to prevent header duplication
    if os.path.exists(output_file):
        os.remove(output_file)

    for i, read in enumerate(bamfile.fetch()):
        if (not read.is_paired or read.is_unmapped or read.mate_is_unmapped 
            or read.mapping_quality < mapq_threshold):
            continue  # Skip low-quality reads

        read_chrom = read.reference_name
        read_pos = read.reference_start
        mate_chrom = read.next_reference_name
        mate_pos = read.next_reference_start

        # **Check for large deletions in CIGAR string**
        read_cigar = read.cigarstring
        deletion_sizes = [int(d) for d in re.findall(r"(\d+)D", read_cigar)] if read_cigar else []
        max_deletion = max(deletion_sizes) if deletion_sizes else 0  # Largest deletion found

        # **If mate is already seen, process the pair**
        if read.query_name in read_pairs:
            mate = read_pairs.pop(read.query_name)

            # **Check mate mapping quality**
            if mate.mapping_quality < mapq_threshold:
                continue  # Skip pairs where mate has low mapping quality

            # **Assign AS score if missing for both read and mate**
            read_alignment_score = read.get_tag("AS") if read.has_tag("AS") else 0
            mate_alignment_score = mate.get_tag("AS") if mate.has_tag("AS") else 0

            # **Collect SA tags for both read and mate**
            sa_tag_read = read.get_tag("SA") if read.has_tag("SA") else None
            sa_tag_mate = mate.get_tag("SA") if mate.has_tag("SA") else None

            # **Classify Chimeric Events**
            if read_chrom != mate_chrom:  # Interchromosomal translocation
                event_type = "Translocation"
            elif abs(read_pos - mate_pos) > distance_threshold:  # Large deletion or inversion
                event_type = "Large Deletion or Inversion"
            elif max_deletion >= deletion_threshold:  # Large deletion detected in CIGAR
                event_type = f"Long Deletion {max_deletion}bp"
            else:
                continue  # Skip non-chimeric reads

            # **Store result in list**
            results.append([read.query_name, read_chrom, read_pos, mate_chrom, mate_pos, 
                            read_cigar, mate.cigarstring, read.mapping_quality, 
                            mate.mapping_quality, read_alignment_score, mate_alignment_score, event_type, 
                            sa_tag_read, sa_tag_mate])

        else:
            read_pairs[read.query_name] = read  # Store read for later pairing

        # **Write to disk periodically to avoid memory overload**
        if i > 0 and i % chunk_size == 0:
            write_header = not os.path.exists(output_file)  # Write header only for the first chunk
            write_chunk_to_file(results, output_file, write_header)
            print(f"Processed {i} reads...")

    # **Write remaining results to file**
    if results:
        write_header = not os.path.exists(output_file)
        write_chunk_to_file(results, output_file, write_header)

    print(f"High-confidence chimeric reads saved to {output_file}")

def main():
    bam_location = sys.argv[1]  # Path to BAM file
    sample = sys.argv[2]
    bamfile_path = f"{bam_location}/sorted_{sample}_aligned.bam"
    
    output_file = f"{sample}_high_confidence_chimeric_reads.tsv"
    extract_high_confidence_chimeric_reads(bamfile_path, output_file, chunk_size=10000, distance_threshold=1000, deletion_threshold=50, mapq_threshold=30, alignment_score_threshold=100)

if __name__ == "__main__":
    main()
