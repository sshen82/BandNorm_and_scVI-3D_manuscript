#!/bin/bash

input_folder=""
output_folder=""
python ./data_conversion/interaction_to_sparse_matrix.py ${input_folder} ${output_folder} "human" 1000000 10

/s/R-3.6.1/bin/Rscript ./topic_modeling/runcisTopic_sparse_CGS.R
