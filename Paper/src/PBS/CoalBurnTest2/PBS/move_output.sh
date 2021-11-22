#!/bin/bash -l

# Combine output into a single file
cd /scratch/ht96/nb9894/slim_SeptGadiTest/
cat ./out_stabsel_pos* >> ./out_stabsel_pos.csv
cat ./out_stabsel_burnin* >> ./out_stabsel_burnin.csv
cat ./out_stabsel_means* >> ./out_stabsel_means.csv
cat ./out_stabsel_opt* >> ./out_stabsel_opt.csv
cat ./out_stabsel_muts* >> ./out_stabsel_muts.csv
cat ./out_stabsel_dict* >> ./out_stabsel_dict.csv

# Delete loose files with seed and model indices: keep the LD matrices
find -regex "./out_stabsel_[a-z]*[0-9]*_*[0-9]*.csv+" -delete
