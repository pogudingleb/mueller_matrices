Runtimes averaged over 20 runs
Times reported in microseconds in the form `mean` ± `stdev`
| Method \ Threads | 1 | 2 | 3 | 4|
|----|----|----|----|----|
| Charpoly avoiding H |44139.8 ± 2912.0 | 30790.5 ± 5417.9 | 23347.0 ± 1687.7 | 25186.9 ± 3503.9|
| Sylvester criterion avoiding H |62065.8 ± 2515.6 | 40685.5 ± 6744.2 | 32218.0 ± 2673.9 | 29132.5 ± 3907.2|
| Eigen's self-adjoint eigensolver |1012055.8 ± 55713.0 | 584151.2 ± 70270.5 | 506391.1 ± 39869.8 | 461919.5 ± 60743.1|
| Eigen's Choletsky |87937.1 ± 5245.3 | 56501.0 ± 10423.4 | 42854.4 ± 1566.0 | 42142.1 ± 5502.6|
| Formula via Pauli matrices |24244.7 ± 2258.5 | 16598.1 ± 3680.4 | 10830.0 ± 478.2 | 12453.5 ± 2700.0|

Versions:
 * Eigen: 3.4.0
 * Compiler: Clang 14.0.0
 * OS: MacOS 12.6

 Machine: 1,6 GHz Dual-Core Intel Core i5, 16 GB RAM
