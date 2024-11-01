# signet_v1
SigNet (Significance Networks): Prioritization of causal genes from GWAS using within-locus and between-locus information

Authors: Zeinab Mousavi, Joel Bader, Johns Hopkins University

Contact: joel.bader@jhu.edu

## Installation overview

We like to have a top-level directory with sub-directories for the github repo, for external data (GWAS data, colocalization data, interaction data, genome reference), and for results. When installation is complete, the directory structure should be similar to this, created using `tree`:

```
signet
├── DATA
│   ├── ENSEMBL
│   │   └── ensembl_gene.txt
│   ├── GWAS
│   │   ├── HR_GWAS_rows.txt
│   │   ├── HR_coloc_genes.txt
│   │   ├── HR_list.txt
│   │   ├── HR_rsid_GRCh38p13.txt
│   │   ├── JT_GWAS_rows.txt
│   │   ├── JT_rsid_GRCh38p13.txt
│   │   ├── PR_GWAS_rows.txt
│   │   ├── PR_coloc_genes.txt
│   │   ├── PR_rsid_GRCh38p13.txt
│   │   ├── QRS_GWAS_rows.txt
│   │   ├── QRS_coloc_genes.txt
│   │   ├── QRS_rsid_GRCh38p13.txt
│   │   ├── QT_GWAS_rows.txt
│   │   ├── QT_coloc_genes.txt
│   │   ├── QT_rsid_GRCh38p13.txt
│   │   ├── exome_genes.txt
│   │   ├── gwas_catalog.txt
│   │   └── omim_genes.txt
│   ├── PROTEIN
│   │   └── human_annotated_PPIs.txt
│   └── TRRUST
│       ├── trrust.human.txt
│       ├── trrust.txt
│       ├── trrust_rawdata.human.tsv
│       ├── trrust_rawdata.human.tsv.1
│       └── trust.txt
├── RESULTS
│   └── GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess
│       └── summary_genestats_manuscript.txt
└── signet_v1
    ├── LICENSE
    ├── README.md
    ├── SAMPLEDATA
    │   ├── HR_coloc_genes.txt
    │   ├── HR_rsid_GRCh38p13.txt
    │   ├── JT_rsid_GRCh38p13.txt
    │   ├── PR_coloc_genes.txt
    │   ├── PR_rsid_GRCh38p13.txt
    │   ├── QRS_coloc_genes.txt
    │   ├── QRS_rsid_GRCh38p13.txt
    │   ├── QT_coloc_genes.txt
    │   ├── QT_rsid_GRCh38p13.txt
    │   ├── exome_genes.txt
    │   └── omim_genes.txt
    ├── counts_postprocessing_manuscript.py
    ├── dist_hist.py
    ├── plot_networks.py
    ├── pval_hist.py
    ├── run_postprocessing.sh
    ├── signet.py
    ├── signet_config.yaml
    ├── signet_v1_conda.yaml
    ├── tables_clipboard.py
    ├── utils_data.py
    ├── utils_plot.py
    └── utils_ppi.py
```

## Installation process

Somewhere in your path, create a top-level directory called `signet` and then
```
cd signet
git clone https://github.com/joelbaderlab/signet_v1.git
mkdir DATA
mkdir RESULTS
```

Then using EKG GWAS data from TopMed as an example, we go do the `DATA` directory and create sub-directories:
```
cd DATA
mkdir ENSEMBL
mkdir GWAS
mkdir PROTEIN
mkdir TRRUST
```

## Data preparation: Ensembl

The Ensembl data provides gene and SNP locations. We used Ensembl biomart through the web interface to select the following features:
```
Gene stable ID
Gene stable ID version
Transcript stable ID
Transcript stable ID version
Gene type
Gene name
Gene start (bp)
Gene end (bp)
Transcript start (bp)
Transcript end (bp)
Transcription start site (TSS)
Chromosome/scaffold name
Gene Synonym
Gene description
```

Then we generated the URL and used `wget` as recommended by Ensembl (see http://useast.ensembl.org/info/data/biomart/biomart_restful.html#wget). The XML below changes the default to include the header and provide a completion stamp. The file we downloaded was about 128 MB.
```
cd ENSEMBL
wget -O ensembl_gene.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1"><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_gene_id_version" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ensembl_transcript_id_version" /><Attribute name = "gene_biotype" /><Attribute name = "external_gene_name" /><Attribute name = "start_position" /><Attribute name = "end_position" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "transcription_start_site" /><Attribute name = "chromosome_name" /><Attribute name = "external_synonym" /><Attribute name = "description" /></Dataset></Query>'
cd ..
```

## Data preparation: GWAS

Next is alphabetical order is the GWAS directory. Each phenotype has three files created, named `PHENO_GWAS_rows.txt`, `PHENO_coloc_genes.txt`, `PHENO_rsid_GRCh38p13.txt`, where we keep the genome version in the name of the SNP locations. We used the most recent meta-analysis as sources for the GWAS data.
```
cd GWAS
wget -O gwas_catalog.txt https://www.ebi.ac.uk/gwas/api/search/downloads/full
head -1 gwas_catalog.txt > HR_GWAS_rows.txt
cat gwas_catalog.txt | awk -F'\t' '$2 == "27798624"' >> HR_GWAS_rows.txt
head -1 gwas_catalog.txt > JT_GWAS_rows.txt
cat gwas_catalog.txt | awk -F'\t' '$2 == "27958378"' >> JT_GWAS_rows.txt
head -1 gwas_catalog.txt > PR_GWAS_rows.txt
cat gwas_catalog.txt | awk -F'\t' '$2 == "30046033"' >> PR_GWAS_rows.txt
head -1 gwas_catalog.txt > QRS_GWAS_rows.txt
cat gwas_catalog.txt | awk -F'\t' '$2 == "27659466"' >> QRS_GWAS_rows.txt
head -1 gwas_catalog.txt > QT_GWAS_rows.txt
cat gwas_catalog.txt | awk -F'\t' '$2 == "24952745"' >> QT_GWAS_rows.txt
cd ..
```

We mapped the positions of the SNPs in the GWAS catalog to the same assembly we used for the gene coordinates. We used the Biomart interface using the lists of SNPs from the GWAS catalog rows as a filter. Pre-generated results are available in `SAMPLEDATA`:
```
cd GWAS
cp ../../signet_v1/SAMPLEDATA/*_rsid_GRCh38p13.txt .
cd ..
```

The colocalization data sets were generated by running `coloc`. For a new phenotype, other colocalization results could be used. The results we generated are available in the `SAMPLEDATA` subdirectory and can be copied over for testing:
```
cd GWAS
cp ../../signet_v1/SAMPLEDATA/*_coloc_genes.txt .
cd ..
```

Genes with other special functional features are listed in additional files. Files with genes with significant assocations with exome variation and with Mendelian evidence are available in `SAMPLEDATA`. These types of lists would have to be generated for a new phenotype.
```
cd GWAS
cp ../../signet_v1/SAMPLEDATA/exome_genes.txt .
cp ../../signet_v1/SAMPLEDATA/omim_genes.txt .
cd ..
```

## Data preparation: Protein-protein interactions

Protein-protein interactions are from the IID database:
```
cd PROTEIN
wget 'http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz'
gunzip human_annotated_PPIs.txt.gz
cd ..
```

## Data preparation: Gene-regulatory interactions

Gene-regulatory interactions are from TRRUST. The TRRUST download does not have headers, so we add them.
```
cd TRRUST
wget 'https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv'
echo -e 'TranscriptionFactor\tTarget\tSign\tPMIDs' > trrust.human.txt
cat trrust_rawdata.human.tsv >> trrust.human.txt
cd ..
```

## Software dependencies

It can be helpful to create a clean environment to run software. We tend to use conda. First, create a new environment and switch to it:
```
conda create --name signet_v1
conda activate signet_v1
```
Then we install dependencies through `conda`. Note that the `yaml` dependency is `pyyaml`, which in turn installs `yaml`.
```
conda install python
conda install pyyaml
conda install numpy
conda install scipy
conda install pandas
conda install matplotlib
conda install python-graphviz
```
The resulting full list of conda packackes is provided as `signet_v1_conda.yaml`. We don't recommend creating the conda environment with this list, however. We usually have better luck just specifying the main dependencies.

## Running SigNet

Everything should run fine from here. Assuming you are in the top-level directory with subdirectories `signet_v1` with the cloned repository, `DATA` with the data, and `RESULTS` for results:
```
conda activate signet_v1
cd signet_v1
python signet.py
```

Underneath the `RESULTS` directory will be a subdirectory named `GWAS_PHENO_250kb_initbestguess` where PHENO is a concatenated list of the phenotypes being analyzed together, `250kb` is the flank length, and `initbestguess` indicates that configurations start with the best guess gene as the active gene, rather than an initial random gene selection.

The repository is set to run 5 independent restarts. Although the best guess initialization should be identical for each restart, changing the order of traversing the loci can change the results from restart to resart. The parameter `-s1` sets the number of restarts. For example,
```
python signet.py -s1 10
```
should perform 10 random restarts. The restart number is used as a seed using `random.seed(seed)` in `run_experiment` in `signet.py`. This can help in testing for consistency. Changing to `random.seed()` would use the system time for the random seed.

Running the first time for a particular set of GWAS loci will take longer because the full list of interactions must be scanned. The software extracts the subnet corresponding to the GWAS loci for subsequent runs.

## Summarizing results

Separate programs read the output from the `RESULTS` directory and generate summary tables and plots.

The summary table `summary_genestats_manuscript.txt` is the main file that aggregates the final gene selections across the independent restarts. It is generated by `counts_postprocessing_manuscript.py`:
```
python ./counts_postprocessing_manuscript.py --results_dir ../RESULTS/GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess --seed_start 0 --seed_end 10
```

Histograms of gene weights and selection probabilities are generated by `pval_hist.py`:
```
python pval_hist.py --results_dir ../RESULTS/GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess
```

Histograms of signed distances from the SNP to the selected gene TSS are generated by `dist_hist.py`:
```
python dist_hist.py --results_dir  ../RESULTS/GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess
```

Network views are created by `plot_networks.py`:
```
python plot_networks.py --results_dir  ../RESULTS/GWAS_QT+HR+QRS+PR+JT_250kb_initbestguess
```

These are collected into the shell script `run_postprocessing.sh`.

An additional file, `tables_clipboard.py`, has a clipboard of code that was used previously to generate tables. This may be helpful to see how additional numbers can be generated. Paths in this file are hard-coded and would have to be set properly using the logic of the scripts called by `run_postprocessing.sh` as a template.

## Pathway enrichment

We used Enrichr from the Ma'ayan lab for pathway enrichment, available at https://maayanlab.cloud/Enrichr/

Lists of genes to submit can be generated from the `summary_genestats_manuscript.txt` file:
```
awk -F'\t' '$21 == "signet" {print $18}' < summary_genestats_manuscript.txt > genes_signet.txt
```
Our analysis used `Pathways` -> `KEGG 2021 Human`.

