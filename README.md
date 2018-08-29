# sampling-benchmarks

This repository consists of weighted and unweighted sampling benchmarks. 

The format of these files follow DIMACS format. We extend DIMACS format to allow specification of "Sampling Set" and "Weights" as follows. 

Every benchmark has few lines with "c ind" prefix to specify sampling set. 
For example if the sampling set is {10, 13, 15, 16, 25, 28, 39, 41, 43, 45, 5, 53, 6, 69, 78, 9, 93}

```c ind 10 13 15 16 25 28 39 41 43 45 0  ```</br>
```c ind 5 53 6 69 78 9 93 0```

Note that each "c ind" line lists only 10 variables at most, followed by a closing '0'.

Literal weights are specified  via lines such as the following:

```w 23 0.3```</br>
```w -23 0.4```

These indicate that the positive literal of variable 23 has weight 0.3, and the
negative literal has weight 0.4. Any positive real numbers may be given as weights.
