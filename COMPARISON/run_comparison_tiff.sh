#!/usr/bin/env bash

resultsdir=.
summaryfile=${resultsdir}/s1table_summary_genestats_manuscript_pops.txt
python comparison_signet_pops_tiff.py --resultsdir ${resultsdir} --summaryfile ${summaryfile}
