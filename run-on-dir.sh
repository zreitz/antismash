#!/bin/bash

for file in ../test-genomes/reference/*.gbff; do
  /Users/zach/miniconda3/envs/antismash/bin/python run_antismash.py $file --output-dir ../reference/$(basename $file .gbff) -v --minimal --hmmdetection-strictness strict
done

