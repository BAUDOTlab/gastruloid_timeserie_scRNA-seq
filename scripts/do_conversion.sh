#!/bin/bash

CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate convertSO_to_h5ad


for i in $(ls ../../seuratAnalysis/gastruloidsAnalysis/rdsObjects/06_clusters_lab_day_*)
do
	filename=$(basename $i .rds)
	echo $filename
	Rscript ./convertSO_to_h5ad.R $i ../../lineageInference/inputData/h5adObjects/$filename.h5Seurat
done

conda deactivate
