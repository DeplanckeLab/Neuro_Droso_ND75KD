# Start from docker image from Aerts' lab
FROM aertslab/pyscenic:0.12.1

# Update
RUN apt update && apt upgrade -y

# Upgrade pip and install anndata and ipykernel for making a Jupyter kernel
RUN python -m pip install --upgrade pip anndata==0.11.4 ipykernel==6.10.0 ipyevents==2.0.4 ipywidgets==8.1.8 jupyter-client==7.3.4
