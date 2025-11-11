#!/bin/bash
# Zip 100 directories at a time into batch_###.zip files

nzip=${1:-100}


count=0
batch=1
dirs_to_zip=()

# Sort directories numerically
for d in $(ls -d [0-9][0-9][0-9][0-9] | sort -n); do
    dirs_to_zip+=("$d")
    ((count++))

    # When 100 directories collected, zip them
    if (( count == "$nzip" )); then
        zip_name=$(printf "batch_%03d.zip" "$batch")
        echo "Creating $zip_name ..."
        zip -r -q "$zip_name" "${dirs_to_zip[@]}"
        ((batch++))
        count=0
        dirs_to_zip=()
    fi
done

# Handle remaining directories (less than 100)
if (( count > 0 )); then
    zip_name=$(printf "batch_%03d.zip" "$batch")
    echo "Creating $zip_name ..."
    zip -r -q "$zip_name" "${dirs_to_zip[@]}"
fi

echo "✅ Done: created zips by 100 directories each."

