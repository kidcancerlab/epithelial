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

mkdir -p input/GSE151374
cd input/GSE151374/

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE151374&format=file" \
    -O GSE151374_RAW.tar

tar -xvf GSE151374_RAW.tar
rm GSE151374_RAW.tar
