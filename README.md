# Description
This repository contains the code to recreate the simulation study in Baybutt and Navjeevan (2023) "Doubly-Robust Inference for Conditional Average Treatment Effects with High Dimensional Controls."  This code can be used to replicate the simulation results, but we also hope that it can be useful for practioners who are interested in implementing the doubly-robust inference procedure described in the paper. For those interested in the latter, the main function of interest is `gHatBN(data)` which takes in "data", a list containing the following:

- "Y": A $n \times 1$ vector containing the outcome variable
- "D": A $n \times 1$ vector containing the treatment indicator
- "Z": A $n \times p$ vector containing the high-dimensional covariates needed for conditional ignorability. The first column of "Z" should be a vector of ones.
- "Pmat": A $n \times K$ vector containing the second-stage basis terms.

The function `gHatBN(data)` takes in this list and returns the following objects of interest

- "beta": The $K \times 1$ vector $\hat\beta$
- "gHat": A vector of predicted conditional counterfactual outcomes obtained by multiplying "Pmat" by "beta"
- "Omega": The matrix \Omega from Section 4 of the paper. This matrix governs the asymptotic variance of the final series estimator and can be used for inference.
