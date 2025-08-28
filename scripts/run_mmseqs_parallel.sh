#!/bin/bash

# SK_DB clustering
mmseqs cluster SK_DB SK_DB_clu_id75 SKtmp_id75 --min-seq-id 0.75 && \
mmseqs createtsv SK_DB SK_DB SK_DB_clu_id75 SK_DB_clu_id75.tsv &

mmseqs cluster SK_DB SK_DB_clu_id90 SKtmp_id90 --min-seq-id 0.9 && \
mmseqs createtsv SK_DB SK_DB SK_DB_clu_id90 SK_DB_clu_id90.tsv &

# CL_DB clustering
mmseqs cluster CL_DB CL_DB_clu_id75 CLtmp_id75 --min-seq-id 0.75 && \
mmseqs createtsv CL_DB CL_DB CL_DB_clu_id75 CL_DB_clu_id75.tsv &

mmseqs cluster CL_DB CL_DB_clu_id90 CLtmp_id90 --min-seq-id 0.9 && \
mmseqs createtsv CL_DB CL_DB CL_DB_clu_id90 CL_DB_clu_id90.tsv &

# CL_SK_DB clustering
mmseqs cluster CL_SK_DB CL_SK_DB_clu_id50 CL_SKtmp_id50 --min-seq-id 0.5 && \
mmseqs createtsv CL_SK_DB CL_SK_DB CL_SK_DB_clu_id50 CL_SK_DB_clu_id50.tsv &

mmseqs cluster CL_SK_DB CL_SK_DB_clu_id75 CL_SKtmp_id75 --min-seq-id 0.75 && \
mmseqs createtsv CL_SK_DB CL_SK_DB CL_SK_DB_clu_id75 CL_SK_DB_clu_id75.tsv &

mmseqs cluster CL_SK_DB CL_SK_DB_clu_id90 CL_SKtmp_id90 --min-seq-id 0.9 && \
mmseqs createtsv CL_SK_DB CL_SK_DB CL_SK_DB_clu_id90 CL_SK_DB_clu_id90.tsv &

# Wait for all background jobs to finish
wait

echo "All clustering and TSV generation completed!"
