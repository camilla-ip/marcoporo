# Test a subset of a 2D run
python -m pdb \
./marcoporo.py extractone \
-bin /well/bsg/microbial/marc/phase2/marcoporo/v1.0 \
-profile None \
-config /well/bsg/microbial/marc/phase2/marcp2/data/00-config/config.txt \
-exptid P2-Lab6-R1-2D \
-indir /well/bsg/microbial/marc/phase2/marcp2/data/01-fast5/P2-Lab6-R1-2D \
-instanceN 000 \
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs True \
-fastq False \
-model False \
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
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs True \
-fastq True \
-model False \
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
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs False \
-fastq False \
-model False \
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
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs True \
-fastq False \
-model False \
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
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs True \
-fastq True \
-model False \
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
-outdir /well/bsg/microbial/marc/phase2/marcp2/data/02-extract \
-pairs True \
-fastq True \
-model False \
-stats True \
-samplesize 10 \
-fastqheaderformat concise

