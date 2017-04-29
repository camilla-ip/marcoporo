#!/usr/bin/env bash

./marcoporo.py runmeta \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-experiments /well/bsg/microbial/marc/phase2/marcp2/data/00-config/experiments.txt \
-samplesize 100 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-exptconstants

./marcoporo.py runmeta \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-experiments /well/bsg/microbial/marc/phase2/marcp2/data/00-config/experimentsR9.txt \
-samplesize 100 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-exptconstants

python -m pdb \
./marcoporo.py exptconstants \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-experiments /well/bsg/microbial/marc/phase2/marcp2/data/00-config/experimentsR9.txt \
-samplesize 10 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract

python -m pdb \
./marcoporo.py extract \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-experiments /well/bsg/microbial/marc/phase2/marcp2/data/00-config/experimentsR9.txt \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq False \
-model False \
-stats False
