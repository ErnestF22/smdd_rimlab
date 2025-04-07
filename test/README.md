- create_corridor_log.cpp
Tests creation of a benchmark-sample corridor log

- create_corridor_log_moving.cpp
Tests creation of a benchmark-sample moving corridor log (simulates motion of an Vehicle along a perfectly sampled corridor)

- estimate_pair_hessian_ise_rtc.cpp
Test multiple scan registration based on ISE (TODO)

- plot_eigs.cpp
Plot eigenvalues and eigenvectors of FIM (Eigen library-based computation)

- plot_eigs_gsl.cpp
Plot eigenvalues and eigenvectors of FIM (EGSL-based computation)

- test_create_corridor_log.cpp
Tests FIM of a benchmark sample created corridor log

- test_fisher_matrix_2d_vs_3d.cpp
Test that compares 2D formulation of FIM and its 3D formulation when applied on a 2D scan

- test_fisher_matrix_3d.cpp
Attempt of extending FIM to 3D case

- test_fisher_matrix_gsl_eigen.cpp
Test that compares eigen-based computation of FIM and its EGSL equivalent

- test_fisher_matrix_rimlab
Test that runs an Eigen-based version of Censi's FIM computation, using an outer-product-based method for computing covariance matrices and related normals

- test_fisher_matrix
Test that runs Censi's FIM computation with original support from GSL/EGSL libraries -> DEPRECATED in Ubuntu > 20.4

- test_gme.cpp
Test Gaussian Mixture Estimator

- test_ise_registration.cpp
Test multiple scan registration based on ISE (TODO)

- test_von_mises_mixture.cpp
Test VMM on CSV input

