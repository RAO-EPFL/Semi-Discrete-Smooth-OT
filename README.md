# Semi-Discrete-Smooth-OT
Semi-Discrete Optimal Transport: Hardness, Regularization and Numerical Solution

Authors: Bahar Taşkesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn

## Contents
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Numerical Experiment](#numerical-experiment)


## Introduction
Semi-discrete optimal transport problems, which evaluate the Wasserstein distance between a discreteand a generic (possibly non-discrete) probability measure, are believed to be computationally hard. Even though such problems are ubiquitous in statistics, machine learning and computer vision, however, this perception has not yet received a theoretical justification. To fill this gap, we prove that computing the Wasserstein distance between a discrete probablity measure supported on two points and the Lebesgue measure on the standard hypercube is already #P-hard. This insight prompts us to seek approximate solutions for semi-discrete optimal transport problems. We thus perturb the underlying transportation cost with an additive disturbance governed by an ambiguous probability distribution, and we introduce a distributionally robust dual optimal transport problem whose objective function is smoothed with the most adverse disturbance distributions from within a given ambiguity set. We further show that smoothing the dual objective function is equivalent to regularizing the primal objective function, and we identify several ambiguity sets that give rise to several old and new regularization schemes. As a byproduct, we discover an intimate relation between semi-discrete optimal transport problems and discrete choice models traditionally studied in psychology and economics. To solve the regularized optimal transport problems efficiently, we use a stochastic gradient descent algorithm with imprecise stochastic gradient oracles. A new convergence analysis reveals that this algorithm improves the best known convergence guarantee for popular semi-discrete optimal transport problems with entropic regularizers.

## Quick Start
This repository contains of computation of smooth semi-discrete optimal transport problem using Algorithm 1 that is presented in the paper. Install YALMIP from https://yalmip.github.io/tutorial/installation/ and MOSEK from https://docs.mosek.com/9.2/install/installation.html by following the instructions.


## Numerical Experiment
The results in the numerical experiments section (Figure 1) are obtained by running main.m script.
All the supplementary functions are placed under the src folder.




