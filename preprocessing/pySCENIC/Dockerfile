# Start from docker image from Aerts' lab
FROM aertslab/pyscenic:0.12.1

# Update
RUN apt update && apt upgrade -y

# Upgrade pip and install anndata (h5ad loading), ipykernel (jupyter), regdiffusion (to replace grnboost2)
RUN python -m pip install --upgrade pip anndata==0.10.9 pandas==1.5.1 numpy==1.23.5 regdiffusion==0.1.1 asttokens==2.4.1 comm==0.2.1 debugpy==1.8.1 decorator==5.1.1 et-xmlfile==1.1.0 exceptiongroup==1.2.0 executing==2.0.1 jedi==0.19.1 matplotlib-inline==0.1.6 openpyxl==3.1.2 packaging==21.3 parso==0.8.3 platformdirs==4.2.0 prompt-toolkit==3.0.43 pure-eval==0.2.2 pygments==2.17.2 pyzmq==25.1.2 tornado==6.1 ipython==8.22.2 traitlets==5.14.1 typing_extensions==4.4.0 wcwidth==0.2.13 ipykernel==6.10.0 ipyevents==2.0.4 ipywidgets==8.1.8 jupyter-client==7.3.4
