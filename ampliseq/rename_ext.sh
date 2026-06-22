for file in *.fastq.gz; do
    # Check if files matching the pattern actually exist
    [ -e "$file" ] || continue

    # 2. Extract the two-part extension (e.g., fastq.gz or fq.gz)
    ext1="${file##*.}"          # gz
    without_ext="${file%.*}"   # strips .gz
    ext2="${without_ext##*.}"   # fastq or fq
    
    # Combine them to get the exact extension format
    # and add the _001 required for ampliseq
    
    #new_ext="_001.${ext1}.${ext2}"   # fastq.gz
    #full_ext="${ext1}.${ext2}"   # fastq.gz
    full_ext="fastq.gz"
    new_ext="_001.fastq.gz"
    
    # 3. get basename
    base="${file%.fastq.gz}"
    

    # 4. Generate the final target filenames
    new_file="${base}$new_ext"
    echo $new_file
    
    # execute mv
    mv -- "$file" "$new_file"
done
