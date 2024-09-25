# Neuro_Droso_ND75-KD
Material and source code for the Neuro_Droso_ND75-KD manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-XXXXX](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXXX). It comprises two 10x libraries (Pan_neuro_control and Pan_neuro_ND75KD).<br/>

## 1. Preprocessing of the single-nuclei libraries

### 1.1. Running CellRanger on the 10x library
Here we follow standard protocols. 

```bash
cellranger count --id Pan_neuro_control --fastqs=./fastq --sample=Pan_neuro_control --transcriptome=${10x_genome} --nosecondary
cellranger count --id Pan_neuro_ND75KD --fastqs=./fastq --sample=Pan_neuro_ND75KD --transcriptome=${10x_genome} --nosecondary
```
**Note:** In the manuscript, we used the dm6 genome assembly.

### 1.2. Generating individual Seurat object for each library <sub>(see full code here: Pan_neuro_control [[Rmd](preprocessing/S0a_Pan_neuro_control.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/Neuro_Droso_ND75KD/blob/main/preprocessing/S0a_Pan_neuro_control.html)] and Pan_neuro_ND75KD [[Rmd](preprocessing/S0b_Pan_neuro_ND75KD.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/Neuro_Droso_ND75KD/blob/main/preprocessing/S0b_Pan_neuro_ND75KD.html)])</sub>

```bash
Rscript -e "rmarkdown::render('./preprocessing/S0a_Pan_neuro_control.Rmd', output_file = './preprocessing/S0a_Pan_neuro_control.html')"
Rscript -e "rmarkdown::render('./preprocessing/S0b_Pan_neuro_ND75KD.Rmd', output_file = './preprocessing/S0b_Pan_neuro_ND75KD.html')"
```

### 1.4. Filtering outlier cells <sub>(see full code here: [[Rmd](code/1.4_Filtering_outlier_cells.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.4_Filtering_outlier_cells.html)])</sub>
At step 1.3 we create a Seurat object containing the raw data counts and a cell metadata containing their assigned TFs. At this stage, following the previous pipeline, all cells have an assigned TF.<br/>
Now, following the standard Seurat pipeline, we aimed at removing outlier cells given the following criteria:
- Remove outlier cells using the `isOutlier` function of the `scater` package for library depth (nCounts) and number of detected genes (nFeature)
- Remove cells with too much mitochondrial RNA or ribosomal RNA
- Remove cells with not enough protein-coding RNA

Here is the function that is performing this filtering:

```R
filtering_outlierCells <- function(seurat_path, libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75){
data.seurat <- readRDS(seurat_path) # Read Seurat object generated at step 1.3

### ---Features and Library size
# Looking for outlier cells in the nCount_RNA distribution (returns TRUE/FALSE array)
libsize.drop <- scater::isOutlier(data.seurat$nCount_RNA, nmads=libsize_nmads, type="lower", log=TRUE) # nCount_RNA / colSums(data.seurat)
# Looking for outlier cells in the nFeature_RNA distribution (returns TRUE/FALSE array)
features.drop <- scater::isOutlier(data.seurat$nFeature_RNA, nmads=features_nmads, type="lower", log=TRUE) # nFeature_RNA / as.vector(colSums(data.seurat > 0))

### --- Gene annotation
data.gene_annot <- fread("metadata/GRCm38_gene_annot.tsv", data.table = F)
  
### --- Mitochondrial
mito.genes <- subset(data.gene_annot, is_mitochondrial)$ensembl_id
mito.genes <- mito.genes[mito.genes %in% rownames(data.seurat)]
data.seurat$percent.mito <- data.seurat[mito.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
### --- Ribosomal
ribo.genes <- subset(data.gene_annot, is_ribosomal)$ensembl_id
ribo.genes <- ribo.genes[ribo.genes %in% rownames(data.seurat)]
data.seurat$percent.rRNA <- data.seurat[ribo.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
### --- Protein Coding
protCod.genes <- subset(data.gene_annot, biotype == "protein_coding")$ensembl_id
protCod.genes <- protCod.genes[protCod.genes %in% rownames(data.seurat)]
data.seurat$percent.ProtCod <- data.seurat[protCod.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100

### Running the filtering
data.seurat <- data.seurat[,features.drop < 1 & libsize.drop < 1 & data.seurat$percent.mito < max_pc_mito & data.seurat$percent.rRNA < max_pc_rRNA & data.seurat$percent.ProtCod > min_pc_protCod]
saveRDS(data.seurat, file = "...") # Saving the filtered Seurat object
}
```

We did not use the same filtering thresholds for all experiments (exp05-12), because of different qualities and biological contexts across experiments. Here is a summary of our filterings:
```R
filtering_outlierCells("exp05", libsize_nmads = 4, features_nmads = 4, max_pc_mito = 10, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp06", libsize_nmads = 4, features_nmads = 4, max_pc_mito = 10, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp07", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 25, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp08", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp09", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 35,  min_pc_protCod = 75)
filtering_outlierCells("exp10", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp11", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp12-13", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 30, max_pc_rRNA = 60,  min_pc_protCod = 75)
```

## 2. Manuscript Figures

All manuscript figures are reproducible through scripts deposited [here](./figures/)
