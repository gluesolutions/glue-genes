![glue_genes_icon_big](https://user-images.githubusercontent.com/3639698/135276631-e67b413e-9299-4bb8-be68-fa18f06759b9.png)

The glue-genes meta-package
=========================



This package does not contain any code, but instead is a meta-package that 
installs all the packages necessary to use glue with genomics data. 

Due to some complicated binary dependencies, the recommended procedure for 
installing this package is to install into a new conda environment as follows:

```
conda create -n glue-genes python==3.8
conda activate glue-genes
conda install -c glueviz glueviz
conda install -c bioconda pairix
conda install -c bioconda tabix
pip install glue-genes@git+https://github.com/gluesolutions/glue-genes.git
```
