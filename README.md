# tswHist

**tswHist** (pronounced "tswist" like "twist", but with an additional 's') stands for **Turbo Sliding Window Histograms**.  
It is a lightweight MATLAB + C library for fast, efficient computation of sliding-window histograms on 1D signals.  
The project provides both pure MATLAB and high-performance MEX (C) implementations.

---

## Features

- Efficient sliding window histogram computation for large 1D signals
- Differential update algorithm for speed using [S. Perreault and P. Hebert, "Median Filtering in Constant Time,"](https://doi.org/10.1109/TIP.2007.902329)
- Multiple variants: pure MATLAB, custom MATLAB, and MEX (C) backends
- Test and benchmarking

---

## File Overview

| File/Folder               | Description                                                                                 |
|---------------------------|---------------------------------------------------------------------------------------------|
| `tswHist.m`               | Main MATLAB function for sliding window histograms (multiple algorithm variants)            |
| `tswHist_mx.c`            | Twin MEX function for tswHist.m                                                             |
| `tswHist_mx.h`            | C header with core routines for `tswHist_mx.c` and `hist_int_mx.c`                          |
| `hist_int_mx.c`           | Twin MEX function for local hist_int matlab function (used by `tswHist.m` custom-mx variant)|
| `Makefile`                | Build script for compiling all MEX files                                                    |
| `test/test_tswHist.m`     | Test script for validating correctness and benchmarking all implementations                 |

---

## Requirements

- MATLAB
- GNU make

(tested on Debian linux)

## Build

To build the MEX files (`tswHist_mx`, `hist_int_mx`):

```sh
make
```

For debug builds:
```sh
make debug
```

## Usage
In MATLAB, add the project folder to your path and use:

```matlab
[histMat, loci, edges] = tswHist(x, n_bins, win_len, stride, variant)
```

* `x`: Input vector (1D signal)
* `n_bins`: Number of histogram bins
* `win_len`: Sliding window length
* `stride`: Step size for sliding window
* `variant`: (optional) 'builtin', 'custom-ml', or 'custom-mx' (default: 'builtin')

For maximum speed, use the MEX implementation directly:

```matlab
[histMat, loci, edges] = tswHist_mx(x, n_bins, win_len, stride)
```

## Testing
Run the test script to validate functionality and performance:

```sh
make test
```

On my computer, I get the following results:

| Implementation                                                               | Execution time (s) |
|------------------------------------------------------------------------------|--------------------|
| Exhaustive computation for each window                                       | 0.4423             |
| `tswHist` using builtin `histcounts` for 1st window histogram                | 0.0158             |
| `tswHist` using personal MATLAB histogram for 1st window histogram             | 0.0172             |
| `tswHist` using personal MEX histogram for 1st window histogram                | 0.0166             |
| `tswHist_mx` (C MEX twin of `tswHist.m`)                                     | 0.0055             |

*Note: Lower execution time is (obviously) better.


## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

## References

- [FastHist: https://github.com/cyber-g/FastHist](https://github.com/cyber-g/FastHist)
- [Median Filtering in Constant Time, Perreault & Hebert, IEEE TIP 2007](https://doi.org/10.1109/TIP.2007.902329)