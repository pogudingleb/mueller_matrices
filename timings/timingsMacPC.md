Runtimes averaged over 20 runs
Times reported in microseconds in the form `mean` ± `stdev`
| Method \ Threads | 1 | 2 | 3 | 4 | 5 | 6|
|----|----|----|----|----|----|----|
| Charpoly avoiding H |26828.7 ± 149.6 | 13495.8 ± 97.7 | 9143.4 ± 139.7 | 6851.4 ± 41.0 | 5637.8 ± 44.0 | 4720.8 ± 31.8|
| Sylvester criterion avoiding H |34600.7 ± 263.2 | 17703.3 ± 100.1 | 11905.0 ± 51.4 | 8986.1 ± 49.4 | 7393.9 ± 21.8 | 6259.9 ± 238.3|
| Eigen's self-adjoint eigensolver |561293.2 ± 1654.6 | 289101.2 ± 724.3 | 194285.4 ± 3019.4 | 145544.2 ± 381.7 | 119189.8 ± 555.7 | 99509.6 ± 314.8|
| Eigen's Choletsky |49587.2 ± 435.0 | 25492.0 ± 42.7 | 17073.4 ± 101.6 | 12956.7 ± 48.5 | 10567.1 ± 32.7 | 8827.0 ± 38.1|
| Formula via Pauli matrices |13102.5 ± 75.2 | 6635.9 ± 37.9 | 4532.1 ± 28.0 | 3503.8 ± 41.6 | 2941.1 ± 22.8 | 2568.0 ± 172.8|

Versions:
 * Eigen: 3.4.0
 * Compiler: Clang 13.1.6
 * OS: MacOS 12.3

 Machine: 3.1 GHz 6-Core Intel Core i5, 16 GB RAM
