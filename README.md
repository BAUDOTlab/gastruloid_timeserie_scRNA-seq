# Time series single-cell RNA sequecing data analysis

This repository contains the rmarkdown scripts and python notebooks used to perform the scRNA-seq data analysis presented in Argiro et al.

## Data

### 1. Experimental data

### 2. Atlas reference data

## Seurat analysis

Before running the following command, please go to the seurat directory:

```
cd scripts/seurat
```

### 1. Preprocess the reference atlas of Pijuan-Sala et al.

```
conda activate seuratV4
Rscript -e "rmarkdown::render('atlasExploration.Rmd', output_file = '../../reports/atlasExploration.html')"
```

### 2. Run a standard scRNA-seq analysis on each dataset independently

As the analysis of a single dataset is very similar to the analysis of integrated datasets,
a common .R script has been developped to run either the `singleDataset_pipeline.Rmd`
or the `combinedDatasets_pipeline.Rmd`.

To investigate the parameters, enter the command in a terminal

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript pipeline_launcher.R -h
```

The analyses of each single dataset described by the Figure 2 of the article were executed with the following commands:

```
#conda activate seuratV4		# not necessary if the environment is already activated

# lab_day_04
Rscript pipeline_launcher.R -N modelGastruloid -i lab_day_04 -m 7.5 --min_counts 5000 -v vst -w 2000

# lab_day_05
Rscript pipeline_launcher.R -N modelGastruloid -i lab_day_05 -m 7.5 -q 20 --min_counts 2500 -v vst -w 2000

# lab_day_06
Rscript pipeline_launcher.R -N modelGastruloid -i lab_day_06 -q 15 --min_counts 2500 -v vst -w 2000

# lab_day_11
Rscript pipeline_launcher.R -N modelGastruloid -i lab_day_11 -q 8 --min_counts 2500 -v vst -w 2000
```

### 3. Look for markers between Cardiomycytes clusters of day 11

The focus on the cardiomycytes described in the Figure 3.a-f of the article were
created with the `cardiomyocytesMarkers_day11_cl.6-8-9.Rmd` script.

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript -e "rmarkdown::render('cardiomyocytesMarkers_day11_cl.6-8-9.Rmd', output_file = '../../reports/cardiomyocytesMarkers_day11_lab_day_11_cl.6-8-9.html')"
```

### 4. Look for markers between clusters 16 and 18 of day 11

The previous code also generated the Figure 4.a and 4.c-m of the article. The
Figure 4.n was generated with the following code:

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript -e "rmarkdown::render('cluster16_zoomIn.Rmd', output_file = '../../reports/cluster16_zoomIn.html')"
```

### 5. Run a standard scRNA-seq analysis on the merged datasets

After analyzing each single dataset, we merged them using the following parameters:

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript pipeline_launcher.R -i lab_4_days --min_counts 2500 -t 15 -N ciml5  -v vst -w 2000 --merge
```

This command allowed us to create the Figure 5.a of the article.

## Data preparation to use in other formats

In order to use cellRank, a Python package, we created a dataset starting from
the row data, trimmed from the cells that were removed during quality control and
doublets detection. The following code realized this step:

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript -e "rmarkdown::render('dataPreparation_to_PythonAnalysis.Rmd', params = list(dataset = 'lab_4_days'), output_file = '../../reports/cellrank_dataPreparation_lab_4_days.html')"
```

A new .Rds file is created and called `99_rawData_filtered_lab_4_days.rds`.

## cellRank analysis

### 1. Convert the .Rds file into a .h5ad file

As cellRank is a Python package, the .Rds object obtained from the previous step
needs to be converted into a .h5ad object. This object will be loaded as an 
AnnData object.

Due to some package compatibility issues, we also need to change the environment.

```
conda activate convertSO_to_h5ad
cd ../cellRank
Rscript ./convertSO_to_h5ad.R 99_rawData_filtered_lab_4_days.rds 99_rawData_filtered_lab_4_days.h5Seurat
```

The `convertSO_to_h5ad.R` script will convert the .Rds file into a .h5Seurat file.
The .h5Seurat file will then be converted to a .h5ad file, readable format for Python.


### 2. Run the cellRank analysis

To run the Jupyter notebook, the `veloRossi` environment was set. The Figure 5.b
and Figure 5.c were obtained from the cellRank analysis.

```
conda activate veloRossi
jupyter nbconvert --execute --to html cellRank_WOT_lab_4_days_scpMerge.nbconvert.ipynb
```

## URD analysis

URD is an .R package that uses URD objects. First, we export the seurat object into
an URD object. Then, we ran the analysis.

### 1. Change object format

We load the seurat object from the file `99_rawData_filtered_lab_4_days.rds`.
UMAP coordinates and metadata are added to the seurat object. Then, the seurat
object is converted into an urd object.

```
#conda activate seuratV4		# not necessary if the environment is already activated
Rscript exportToURD.R
```

The file `urdObject_lab_4_days.rds` is created.


### Run the URD analysis

Due to some package compatibility issues, we set the `urd` environment. We obtained
the figures Fig.5d-o of the article.

```
conda activate urd
Rscript -e "rmarkdown::render('11.1_knnOutliers.Rmd', output_file = '../../reports/urd_lab_4_days.html')"
```





