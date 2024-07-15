# signet_v1
SigNet (Significance Networks): Prioritization of causal genes from GWAS using within-locus and between-locus information

Authors: Zeinab Mousavi, Joel Bader, Johns Hopkins University

Contact: joel.bader@jhu.edu

## Installation

We like to have a top-level directory with sub-directories for the github repo, for external data (GWAS data, colocalization data, interaction data, genome reference), and for results. So, somewhere in your path, create a top-level directory called `signet` and then
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
```


