## Before compiling

The code uses `Eigen` library for linear algebra in C++. The library should be downloaded from the [website](https://eigen.tuxfamily.org/index.php?title=Main_Page), and the path to library should be inserted to the `EIGEN_PATH` variable in the `Makefile`.

## Compiling

The project can be compiled by typing `make`. 
This will create two binary files `test` and `matrices`.
For subsequent compilations, consider typing `make clean` first.

## Running

Correctness checks can be run by `./test`.

For running the main code filtering the matrices, type
```bash
./matrices method num_threads
```
Here `method` is the index of the method to be used (a full list will be printed if you type just `./matrices` without any arguments) and `num_threads` is the number of threads to be used for the computations.

In order to collect the statistics across different methods and different thread count, one can use a python script `produce_stats.py` as follows:
```bash
python3 produce_stats.py num_thread num_runs
```
Here `num_threads` is the largest number of threads to try (so the code will run each method with `1, 2, ..., num_threads` threads) and `num_runs` is the number of runs to preform for averaging the runtime. The default values are `4` and `20`, respectively.
 
