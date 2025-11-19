# Neuro Droso ND-75KD - pySCENIC
pySCENIC is used for computing regulons from scRNA-seq data.

In this subfolder are all the scripts needed to run the pySCENIC methodology. In particular, we've created a Docker containing all tools/packages to run all the steps of the pipeline and generate the figures. We also created a kernel file, to use the Docker as a Jupyter Notebook kernel.

# Installation
## Docker image
First, we build a local Docker image using the [Dockerfile](Dockerfile) (here we call the local image "pyscenic:0.12.1"). For this, we use the [aertslab/pyscenic:0.12.1](https://hub.docker.com/r/aertslab/pyscenic/tags) Docker image available on the DockerHub, and we add missing dependencies (such as anndata for reading .h5ad files or ipykernel for JupyterHub compatibility).

```bash
docker build -t pyscenic:0.12.1 .
```

## [Optional] Create a JupyterNotebook kernel using the Docker image
Then, you can copy the provided kernel file: [kernel.json](kernel.json) to the Jupyter notebook kernel folder on your local machine.

e.g. on Linux, it is stored in /usr/local/share/jupyter/kernels/

# Running pySCENIC and adding regulon data to the Seurat object
Now, we will run the pySCENIC pipeline in Python, using the docker image that we created before.

## Running pySCENIC in Python <sub>(see full code here: [[ipynb](../S3a_Pan_neuro_integrated_pySCENIC_pipeline.ipynb)])</sub>
Running pySCENIC without any error was far from being an easy journey. We'd like to point an eventual reader to this GitHub issue that we've created for this purpose: [aertslab/pySCENIC/issues/534](https://github.com/aertslab/pySCENIC/issues/534). It explains the issues we faced and how we resolved them. This may help to understand some filtering/mapping that we do in the ipynb file, which are not from the default pySCENIC tutorial.

It is worth noting that we saved all intermediary files in .tsv files. We also applied the binarization method for each regulon, to get a binary output (instead of a continuous regulon score) of activated vs. non-activated regulon. Both regulon outputs are then stored in separated .tsv files to be further imported into the Seurat object.

## Importing pySCENIC regulons to the Seurat object
We do this after all metadata were generated, in step 1.4 of main pipeline
