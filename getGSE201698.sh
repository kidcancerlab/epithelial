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

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

ml purge

mkdir -p input/GSE201698
cd input/GSE201698/

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE201698&format=file" \
    -O GSE201698_RAW.tar

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE201698&format=file&file=GSE201698%5Fcell%5Fbarcode%5Finformation%2Exlsx" \
    -O GSE201698_cell_barcode_information.xlsx

tar -xvf GSE201698_RAW.tar
rm GSE201698_RAW.tar