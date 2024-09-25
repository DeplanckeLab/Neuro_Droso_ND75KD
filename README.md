# Neuro_Droso_ND75-KD
Material and source code for the Neuro_Droso_ND75-KD manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-XXXXX](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXXX). It comprises two 10x libraries (Pan_neuro_control and Pan_neuro_ND75KD).<br/>

## 1. Preprocessing of the single-nuclei libraries

### 1.1. Running CellRanger on the 10x library
Here we follow standard protocols. 

```bash
cellranger count --id Pan_neuro_control --fastqs=./fastq --sample=Pan_neuro_control --transcriptome=${10x_genome} --expect-cells=10000 --chemistry=auto --include-introns=true
cellranger count --id Pan_neuro_ND75KD --fastqs=./fastq --sample=Pan_neuro_ND75KD --transcriptome=${10x_genome} --expect-cells=10000 --chemistry=auto --include-introns=true
```
**Note:** In the manuscript, we used the dm6 genome assembly.

### 1.2. Generating individual Seurat object for each library <sub>(see code & output here: Pan_neuro_control.[[Rmd](preprocessing/S0a_Pan_neuro_control.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/Neuro_Droso_ND75KD/blob/main/preprocessing/S0a_Pan_neuro_control.html)] and Pan_neuro_ND75KD.[[Rmd](preprocessing/S0b_Pan_neuro_ND75KD.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/Neuro_Droso_ND75KD/blob/main/preprocessing/S0b_Pan_neuro_ND75KD.html)])</sub>

```bash
Rscript -e "rmarkdown::render('./preprocessing/S0a_Pan_neuro_control.Rmd', output_file = './preprocessing/S0a_Pan_neuro_control.html')"
Rscript -e "rmarkdown::render('./preprocessing/S0b_Pan_neuro_ND75KD.Rmd', output_file = './preprocessing/S0b_Pan_neuro_ND75KD.html')"
```

### 1.3. Integration of transcriptomics data with Harmony <sub>(see full code here: [[Rmd](code/1.4_Filtering_outlier_cells.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.4_Filtering_outlier_cells.html)])</sub>

[TODO]

## 2. Manuscript Figures

All manuscript figures are reproducible through scripts deposited [here](./figures/)
