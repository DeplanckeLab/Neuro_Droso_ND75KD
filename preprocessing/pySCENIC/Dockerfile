# Start from docker image from Aerts' lab
FROM aertslab/pyscenic:0.12.1

# Update
RUN apt update && apt upgrade -y

# Upgrade pip and install anndata (h5ad loading), ipykernel (jupyter), regdiffusion (to replace grnboost2)
RUN python -m pip install --upgrade pip anndata==0.11.4 pandas==1.5.1 numpy==1.23.5 regdiffusion==0.1.1 ipykernel==6.10.0 ipyevents==2.0.4 ipywidgets==8.1.8 jupyter-client==7.3.4
