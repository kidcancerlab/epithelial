#!/bin/sh
#SBATCH --account=gdrobertslab
#SBATCH --error=slurmOut/getGeo_%j.txt
#SBATCH --output=slurmOut/getGeo_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name getGeo
#SBATCH --wait
#SBATCH --time=0-01:00:00
#SBATCH --partition=general,himem
#SBATCH --mail-user=matthew.cannon@nationwidechildrens.org

set -e ### stops bash script if line ends with error

echo ${HOSTNAME}

ml purge

# want barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
mkdir -p input/GSE129605
cd input/GSE129605/

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129605&format=file" \
    -O GSE129605_RAW.tar

tar -xvf GSE129605_RAW.tar
rm GSE129605_RAW.tar

for file_name in *_barcodes.tsv.gz
do
    base_name=${file_name%%_barcodes.tsv.gz}
    no_number=${base_name%.[0-9]*}
    echo ${base_name}
    echo ${no_number}
    mkdir ${no_number}
    mv ${base_name}*gz ${no_number}/
    rename ${base_name}_ "" ${no_number}/*
    mv ${no_number}/genes.tsv.gz ${no_number}/features.tsv.gz
done




