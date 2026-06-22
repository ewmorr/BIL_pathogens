
# Loop through all forward read files
for file1 in *R1*.fastq.gz; do
    # Check if files matching the pattern actually exist
    [ -e "$file1" ] || continue

    # 1. Identify the corresponding reverse read file
    file2="${file1/R1/R2}"
    
    # Safety check: Ensure the matching R2 file actually exists
    if [ ! -f "$file2" ]; then
        echo "Warning: Matching pair for $file1 not found. Skipping."
        continue
    fi

    # 2. Extract the two-part extension (e.g., fastq.gz or fq.gz)
    ext2="${file1##*.}"          # gz
    without_ext2="${file1%.*}"   # strips .gz
    ext1="${without_ext2##*.}"   # fastq or fq
    
    # Combine them to get the exact extension format
    full_ext="${ext1}.${ext2}"   # fastq.gz

    # 3. Clean the base names (strip extension and replace remaining dots with underscores)
    base1="${file1%.$full_ext}"
    clean_base1="${base1//./_}"
    
    base2="${file2%.$full_ext}"
    clean_base2="${base2//./_}"

    # 4. Generate the final target filenames with R1 and R2 swapped
    # Note: clean_base1 contains "R1", so we swap it to "R2" for its final name
    new_file1="${clean_base1/R1/R2}.$full_ext"
    new_file2="${clean_base2/R2/R1}.$full_ext"

    # 5. Safe execution using a temporary placeholder to prevent clashing
    tmp_file="TMP_SWAP_$file1"
    
    echo "Swapping and cleaning: $file1 <-> $file2"
    mv -- "$file1" "$tmp_file"
    mv -- "$file2" "$new_file2"
    mv -- "$tmp_file" "$new_file1"
done
