Configuration files for the most relevant executables.

Example of usage is: from build folder, call the following bash commmand 

> $ ./test_vomp -cfg ../cfg/gmm_estim.cfg

If one desires to set a parameter differently from what is in the config file, use the same syntax.
E.g., if one wants to set param "in" to "sample.csv", use

> $ ./test_vomp -cfg ../cfg/gmm_estim.cfg -in sample.csv

- gmm_estim.cfg\
    For test_gme executable

- run_vomp_cuda_full.cfg\
    For run_vomp_cuda_full CUDA executable

- test_vomp.cfg\
    For test_vomp executable
