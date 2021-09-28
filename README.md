
The glue-genes meta-package
=========================


![glue-genes](https://user-images.githubusercontent.com/3639698/135114413-8fde24e2-715a-4c66-a8c3-585508266298.png)

This package does not contain any code, but instead is a meta-package that 
installs all the packages necessary to use glue with genomics data. 

Due to some complicated binary dependencies, the recommended procedure for 
installing this package into a new conda environment is as follows:

```
conda create -n glue-genes python==3.8
conda activate glue-genes
conda install -c glueviz glueviz
conda install -c bioconda pairix
conda install -c bioconda tabix
pip install glue-genes@git+https://github.com/gluesolutions/glue-jax.git
```
