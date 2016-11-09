#!/usr/bin/env bash

# CI 2016-11-09
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py exptconstants \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-experiments /well/bsg/microbial/marc/phase2/marcp2/data/00-config/experiments.txt \
-samplesize 300 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract

# CI 2016-11-09
nohup nice \
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P1b-Lab2-R2-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P1b-Lab2-R2-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise \
&> /well/bsg/microbial/marc/phase2/marcp2/data/02-extract/P1b-Lab2-R2-2D_extractone.log &

nohup nice \
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-1D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-1D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise \
&> /well/bsg/microbial/marc/phase2/marcp2/data/02-extract/P2-Lab6-R1-1D_extractone.log &

nohup nice \
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise \
&> /well/bsg/microbial/marc/phase2/marcp2/data/02-extract/P2-Lab6-R1-2D_extractone.log &

nohup nice \
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab7-R1-1D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab7-R1-1D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise \
&> /well/bsg/microbial/marc/phase2/marcp2/data/02-extract/P2-Lab7-R1-1D_extractone.log &

nohup nice \
/well/bsg/microbial/marc/phase2/marcoporo/v1.0/marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab7-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab7-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise \
&> /well/bsg/microbial/marc/phase2/marcp2/data/02-extract/P2-Lab7-R1-2D_extractone.log &

