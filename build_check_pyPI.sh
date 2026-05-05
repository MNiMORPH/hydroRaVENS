#! /bin/sh

# Clean old files in order to not get "400 Client Error: File already exists."
rm -f dist/*.tar.gz
rm -f dist/*.whl

# Build the new distribution
python3 -m build

# Check the new distribution for a successful build
twine check dist/*
