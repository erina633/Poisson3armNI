# Poisson3armNI: R codes for " New Approaches for Testing Non-inferiority for Three-arm Trials with Poisson Distributed Outcomes"

The codes are provided to generate sample size and power for both Frequentist and Bayesian approaches. The following notations are used in this package:

•	lamE (or l_E), lamR (or l_R),  and lamP (or l_P): The Poisson rate parameter in the arm E, R, and P respectively

•	nP: The sample size in the placebo arm (P)

•	nR: The sample size in the reference arm (R)

•	nE: The sample size in the experimental arm (E)

•	theta: effect retention parameter

•	r.alloc (or alloc): Allocation vector which can be (1:1:1), (2:2:1) or (3:2:1) for nE:nR:nP

•	N: Total sample size

We give brief description of the R files below:

1. freq_power_sample

This function calculates the Frequentist sample size for a given value of theta, allocation, alpha, l_P, l_R, and l_E.

Arguments

•	n_total: Maximum number of n to get the sample size

•	alloc: Allocation vector

•	l_E: lambda of the arm E

•	l_R: lambda of the arm R

•	l_P: lambda of the arm P

•	theta: effect retention parameter

Output: Power, power curves for different allocations, power curves for different theta, and minimum sample size of the arm P satisfying power>=1-beta 


2. type1error_fullybayesian.R

This function calculates the estimated type I error under exact Bayesian approach for different values of theta, calculated sample size, allocation, and lamE.

Arguments

•	n: sample size

•	r.alloc: Allocation vector

•	lamE: lambda of the arm E

•	theta: effect retention parameter

Output: Estimated type I error for exact Bayesian approach
