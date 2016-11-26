# Test a subset of a 2D run
python -m pdb \
./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq False \
-stats False \
-samplesize 10

python -m pdb \
./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq True \
-stats False \
-samplesize 10 \
-fastqheaderformat poretools

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs False \
-fastq False \
-stats False \
-samplesize 10 \
-fastqheaderformat concise

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq False \
-stats False \
-samplesize 10 \
-fastqheaderformat concise

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq True \
-stats False \
-samplesize 10 \
-fastqheaderformat concise

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq True \
-stats True \
-fastqheaderformat concise

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P1b-Lab2-R2-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P1b-Lab2-R2-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract \
-pairs True \
-fastq True \
-stats True \
-samplesize 10 \
-fastqheaderformat concise

# Rerunning this one

./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract-test \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise

# Testing just a few files to find the one causing the segmentation fault
./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract-test/indir \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/03-extract-test \
-pairs False \
-fastq True \
-stats True \
-fastqheaderformat concise

