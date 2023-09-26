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
mkdir -p input/GSE132771
cd input/GSE132771/

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132771&format=file" \
    -O GSE132771_RAW.tar

tar -xvf GSE132771_RAW.tar

rm *NML* *IPF* *SCD*
rm GSE132771_RAW.tar

for file_name in *_barcodes.tsv.gz
do
    base_name=${file_name%%_barcodes.tsv.gz}
    no_gsm=${base_name#GSM*_}
    echo ${base_name}
    mkdir ${no_gsm}
    mv *${base_name}* ${no_gsm}/
    rename ${base_name}_ "" ${no_gsm}/*
    mv ${no_gsm}/genes.tsv.gz ${no_gsm}/features.tsv.gz
done




