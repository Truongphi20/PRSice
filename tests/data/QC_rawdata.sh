#!/usr/bin/bash

set -ue

## 1. QC basedata
### 1.1 QC info score
gunzip -c Height.gwas.txt.gz |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > Height.gz

### 1.2 Remove duplicates
gunzip -c Height.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > Height.nodup.gz

### 1.3 Getting ambigious SNPs
gunzip -c Height.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > Height.QC.gz

### 1.4 Clean up
rm -f  Height.gz Height.nodup.gz

## 2. QC targeted data
### 2.1 unzip
unzip EUR.zip

### 2.2 Standard GWAS QC
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    plink1.9\
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC

### 2.3 Cut off heterozygosity
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink:v1.07dfsg-2-deb_cv1 \
    /usr/lib/debian-med/bin/plink\
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC

### 2.4 Heterozygosity rating
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink:v1.07dfsg-2-deb_cv1 \
    /usr/lib/debian-med/bin/plink\
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC

### 2.5 Get valid sample ids
docker run --rm \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  rocker/r-base:4.5.2 \
  Rscript filter_het.R

### 2.6 Mismatching SNPs
docker run --rm \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  rocker/r-base:4.5.2 \
  Rscript mismatching.R

### 2.7 Sex chromosomes
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink:v1.07dfsg-2-deb_cv1 \
    /usr/lib/debian-med/bin/plink\
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC

# Filter sex
docker run --rm \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  rocker/r-base:4.5.2 \
  Rscript filter_sex.R

### 2.8 Relatedness
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    plink1.9\
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC

### 2.9 Generate final QC'ed target data file
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    plink1.9\
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1

### 2.10 Generate covariate file
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    plink1.9\
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --mind 0.02 \
    --pca 10 \
    --out EUR

docker run --rm \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  rocker/r-base:4.5.2 \
  Rscript covariate.R

### Clean up
rm -f EUR.bed EUR.bim EUR.fam\
      .pversion EUR.QC.log EUR.QC.hh EUR.QC.irem EUR.QC.rel.id \
      EUR.hh EUR.a1 EUR.eigenval EUR.irem EUR.log \
      EUR.mismatch EUR.QC.prune.in EUR.QC.prune.out EUR.QC.sexcheck EUR.QC.snplist \
      EUR.QC.valid EUR.valid.sample

### Create a VCF file to watch
docker run -v $PWD:$PWD --rm -w $PWD biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1 \
    plink1.9\
    --bfile EUR.QC \
    --recode vcf \
    --out EUR.QC