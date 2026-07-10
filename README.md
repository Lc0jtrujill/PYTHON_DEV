# Numerical benchmark workspace

This repository contains reproducible homogeneous-platform benchmarks for the fractional Pouzet manuscript. The benchmark branch uses GNU Octave to execute the official MATLAB sources of `flmm2`, `fhbvm`, and `fhbvm2`, together with Octave implementations of P2--PI2, P3--PI3, and P4--PI4.

The timing protocol fixes all numerical-library thread counts to one, performs one warm-up execution, and reports the median and interquartile range of five wall-clock repetitions.
