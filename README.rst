The glue-jax meta-package
========================

This package does not contain any code, but instead is a meta-package that 
installs all the packages necessary to use glue with genomics data. 

Due to some complicated binary dependencies, the recommended procedure for 
installing this package into a new conda environment is as follows:

conda create -n glue-jax python==3.8
conda activate glue-jax
conda install -c glueviz glueviz
conda install -c bioconda pairix
conda install -c bioconda tabix
pip install glue-jax@git+https://github.com/gluesolutions/glue-jax.git
