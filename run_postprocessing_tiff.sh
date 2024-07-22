#!/usr/bin/env bash

resultsdir=../RESULTS/GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess

#python counts_postprocessing_manuscript.py --results_dir ${resultsdir} --seed_start 0 --seed_end 10
python pval_hist_tiff.py --results_dir ${resultsdir}
python dist_hist_tiff.py --results_dir ${resultsdir}
#python plot_networks.py --results_dir ${resultsdir}
