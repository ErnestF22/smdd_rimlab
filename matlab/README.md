Folder containing MATLAB scripts used during the preliminary phases of the project for quicker prototyping/check of the correctness of various formulas. Note that for some scripts, [MANOPT](https://www.manopt.org/) library is required.

- crt_computation.m\
Computation of Radon Transform (RT) self-correlation

- intervalPow2(v1, v2, res)\
The function operates on discretized number i1 = round(v1/res) and i2 = round(v1/res) where res is the resolution.

- pnebi\
Pnebi is an acronym that stands for (double of) product of negative \
exponential (with absolute value of argument) and BesselI function i.e.,\
>2.0 * $\exp{-||x||} \besseli_n(x)$\
where Bessel function of order k is taken into consideration\
Note that Bessel functions are even.

- pnebi0\
Pnebi but here Bessel function of order 0 is taken into consideration

- pnebi1\
Pnebi but here Bessel function of order 1 is taken into consideration

- rt_corr_fourier_coeffs\
RT Correlation Fourier Coefficients

- rt_corr_max\
Maximum of RT correlation

- rt_correlation\
RT correlation

- rt_gmm_eval\
Evaluation of RT on input Gaussian Mixture Model

- startup\
Disable undesired Matlab respose to pressing some keys e.g. ALT+GR

- test_correlation_translation\
Test effect of translation on RT self-correlation

- test_crt\
More general version of test_correlation_translation

- test_rt_corr_max_manopt\
Test search of RT correlation maximum search through MANOPT library

- test_rt_correlation\
Test computation of RT correlation

- test_rt_gmm_eval\
Test computation/evaluation of RT on input Gaussian Mixture Model
