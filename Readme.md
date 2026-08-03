This repository provides the implementation accompanying the research article

A Ramanujan-Scaled Fast Multipole Method for Solving the Secular Equation

The secular equation is a fundamental problem in numerical linear algebra that arises in symmetric rank-one modifications, divide-and-conquer eigensolvers, least-squares problems, and related applications. Direct evaluation of the secular function requires quadratic computational complexity, making it computationally expensive for large-scale problems.

This work presents a Ramanujan-scaled Fast Multipole Method (FMM) that accelerates the evaluation of the secular function and its derivative while improving numerical stability through carefully constructed Ramanujan-based scaling factors. The proposed approach retains the near-linear complexity of the FMM and is coupled with Newton's method for efficient computation of the roots of the secular equation.

Features
Fast Multipole Method implementation for the kernel 1/(x-y)
Ramanujan-based scaling for numerically stable M2L translations
Efficient evaluation of the secular function and its derivative
Newton iteration for solving the secular equation
Numerical experiments and performance comparisons
Reproducible scripts for the results reported in the paper

Please run the test_fmm3.m script from the folder.
