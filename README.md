### Single Measurement Degeneracy Detection

Dependencies (tested):

- CUDA 12.0
- [ROFL](https://github.com/dlr1516/rofl) library, at least *common* module
- Eigen 3.4
- PCL 1.14

Optional: Censi's EGSL (hard to compile on Ubuntu > 20.04)

Executable scripts for running the methods on dataset sequences are in folder app.

Executable scripts for running the methods on single scans are in folder test.

Compiling should require the very standard pipeline

> $ mkdir build && cd build
> cmake ..
> make

CMakeLists is set up so that both *test* and *app* executables will be found inside build/ folder.

Feel free to contact us directly, or to open an issue if compilation/usage goes wrong at any point.

Before running the executable, please refer to README instructions inside *cfg* folder for additional instructions on how to set the parameters.

If you find this project helpful, or use part of it in your works, please cite:

> @article{lodirizzinifontana2022niars,
>  title={Rotation Estimation Based on Anisotropic Angular Radon Spectrum},
>  author={Rizzini, Dario Lodi and Fontana, Ernesto},
>  journal={IEEE Robotics and Automation Letters},
>  volume={7},
>  number={3},
>  pages={7279--7286},
>  year={2022},
>  publisher={IEEE}
> }

or other works from us.

© 2025 Ernesto Fontana, Dario Lodi Rizzini (RIMLab laboratory at the University of Parma).


