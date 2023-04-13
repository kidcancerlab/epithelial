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

# Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259
# Most of the time, what folks seem to do it rename the cellranger files
# and upload them to GEO. We need to rename everything.
# Normally, once you rename them you could just use the normal approach to read
# in the data, but it looks like they somehow transposed the matrix so we need
# to do things manually.

# High Resolution
# want barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
mkdir -p input/GSE141259/high_res

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FHighResolution%5Fgenes%2Etxt%2Egz" \
    -O input/GSE141259/high_res/features.tsv.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FHighResolution%5Fbarcodes%2Etxt%2Egz" \
    -O input/GSE141259/high_res/barcodes.tsv.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FHighResolution%5Frawcounts%2Emtx%2Egz" \
    -O input/GSE141259/high_res/matrix.mtx.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FHighResolution%5Fcellinfo%2Ecsv%2Egz" \
    -O input/GSE141259/high_res/metadata.csv.gz


# Whole Lung
mkdir -p input/GSE141259/whole_lung

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FWholeLung%5Fbarcodes%2Etxt%2Egz" \
    -O input/GSE141259/whole_lung/barcodes.tsv.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FWholeLung%5Fgenes%2Etxt%2Egz" \
    -O input/GSE141259/whole_lung/features.tsv.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FWholeLung%5Frawcounts%2Emtx%2Egz" \
    -O input/GSE141259/whole_lung/matrix.mtx.gz

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141259&format=file&file=GSE141259%5FWholeLung%5Fcellinfo%2Ecsv%2Egz" \
    -O input/GSE141259/whole_lung/metadata.csv.gz

