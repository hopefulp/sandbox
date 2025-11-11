#!/usr/bin/bash

# Organize HDF5 files into directories for Edison upload

for file in *.h5; do
    # Remove the .h5 extension to get the ID
    id="${file%.h5}"

    # Create directory if it doesn’t exist
    mkdir -p "$id"

    # Move the file into its directory
    mv "$file" "$id/"
done

echo "✅ All .h5 files moved into corresponding directories."

