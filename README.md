# poly-bihc
You may run this project with the following main parts. 

## 1. Precomputing symbolic expressions of $phi_{m,n}^{k}(x)$.

Run `make` with makefile in the directory `./symbolic/`, this will produce polynomial files in `./symbolic/poly`, this include `./poly/poly3d` and `./poly/poly2d` symbolic expressions as binary file for biharmonic basis.

## 2. Computing biharmonic coordinates for specific data case.
This require `CUDA` for precomputing. In directory `./compute`, first run `make` with makefile to compile the `CUDA` program, then run as `./compute [PATH_TO_DATA] [Tgt_Degree] [optional: \"compute_mesh\", \"compute_p2p\" (only)]`. 

## 3. Using GUI to having interfare.
First compile GUI with `CMake` in `./gui/`, run as `./main [PATH_TO_DATA]`, have fun.

[Updating...]