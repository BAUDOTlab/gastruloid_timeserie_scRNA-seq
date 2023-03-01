#!/bin/bash

#CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
#source $CONDA_PATH/etc/profile.d/conda.sh
#conda activate SO2h5ad


for i in $(ls ../../../seuratAnalysis/ciml5/rdsObjects/\*lab_{3,4}_days\*)
do
	filename=$(basename $i .rds)
	echo $filename
	Rscript ./convertSO_to_h5ad.R $i ../../../lineageInference/inputData/h5adObjects/$filename.h5Seurat
done

#Rscript ../../scripts/convertSO_to_h5ad.R /mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis/2022-06-09_completeAnalysis_endIntegration_integ/rdsObjects/07_louvainClusters_lab_4_days.rds lab_4_days_integ_RNAassay.h5Seurat

#conda deactivate
