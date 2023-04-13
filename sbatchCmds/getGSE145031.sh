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

cd input/GSE145031/

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE145031&format=file" \
    -O GSE145031_RAW.tar

tar -xvf GSE145031_RAW.tar

sample_list=(GSM4304609_PBS_AT2_Tomato
             GSM4304610_PBS_AT2_nonTomato
             GSM4304611_Day14_AT2_Tomato
             GSM4304612_Day14_AT2_nonTomato
             GSM4304613_Day28_AT2_Tomato
             GSM4304614_Day28_AT2_nonTomato)

#split the files into sub-directories and rename them so we can pull them into R
for sample_name in ${sample_list[@]}
do
    echo ${sample_name}
    mkdir ${sample_name}
    mv ${sample_name}_barcodes.tsv.gz ${sample_name}/barcodes.tsv.gz
    mv ${sample_name}_gene.tsv.gz ${sample_name}/features.tsv.gz
    mv ${sample_name}.mtx.gz ${sample_name}/matrix.mtx.gz
done

rm GSE145031_RAW.tar
