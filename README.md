# gmmL0 - Sparse Inverse Covariance Estimation Problems  

This site contains data files for the paper "On the Solution of L0-Constrained Sparse Inverse Covariance Estimation Problems" by Dzung Phan and Matt Menickelly   

## Authors

Dzung Phan (IBM Research) and Matt Menickelly (Argonne National Laboratory)

Emails: phandu@us.ibm.com and mmenickelly@anl.gov

## Data-set Description

The synthetic Matlab datasets consist of:

n                : the size of random variables 

true_inv_cov_Mat : the ground truth sparse inverse covariance matrix 

true_nnz         : the ground truth number of nonzeros in the inverse covariance matrix

sample_cov_Mat   : the sample covariance matrix.

The real-world Matlab datasets include the sample covariance. 

## Example of using a dataset in Matlab

>> load Chain500.mat

