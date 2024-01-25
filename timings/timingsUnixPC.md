Runtimes averaged over 20 runs
Times reported in microseconds in the form `mean` ± `stdev`
| Method \ Threads | 1 | 2 | 3 | 4 | 5 | 6|
|----|----|----|----|----|----|----|
| Charpoly avoiding H |24740.0 ± 346.2 | 12710.9 ± 78.6 | 8677.9 ± 43.8 | 6615.7 ± 56.3 | 5398.4 ± 34.4 | 4487.9 ± 50.5|
| Sylvester criterion avoiding H |31395.5 ± 179.7 | 16053.5 ± 75.8 | 10976.0 ± 52.1 | 8305.0 ± 44.8 | 6805.5 ± 64.7 | 5682.4 ± 55.2|
| Eigen's self-adjoint eigensolver |426543.5 ± 2865.1 | 215017.0 ± 538.8 | 146306.1 ± 213.5 | 110276.5 ± 279.4 | 89926.3 ± 380.5 | 74886.6 ± 34.0|
| Eigen's Choletsky |95920.1 ± 454.1 | 48315.6 ± 155.2 | 33258.3 ± 93.4 | 25026.5 ± 126.0 | 20491.0 ± 39.4 | 17106.1 ± 55.0|
| Formula via Pauli matrices |11529.2 ± 110.7 | 6053.6 ± 94.5 | 4171.0 ± 28.2 | 3195.0 ± 49.9 | 2632.1 ± 26.3 | 2260.3 ± 46.5|

Versions:
 * Eigen: 3.4.0
 * Compiler: g++ 11.4.1
 * OS: AlmaLinux 9.3

 Machine: Intel(R) Xeon(R) W-1290P CPU @ 3.70GHz, 20 cores
