# LBMD2Q9

This Matlab code aims at solving Lattice Boltzmann Method (LBM) in a 2D channel flow. It is modified from https://github.com/siramirsaman/LBM, with several times of speedup and well-organized codes.

## Goal

This modified version is the first step for myself to understand the basics of LBM.

The goal is to observe how well different LBM methods for capturing a curved wall bounadry condition behave.
The benchmark is chosen to be the drag coefficient obtained from flow over a cylinder, while different methods are applied to capture the curved boundary.

## Performance

I found that in many places vectorization for the innermost loop actually degraded the performance! This is new JIT era.
It is possible, but maybe not worth it, to avoid if statements inside nested loops by using logical indexes. It is relatively easy to implement in Julia.

If the vectorized version is not along the continuous dimension, it makes
much less sense to do so. In fact, in some cases it may even degrade the performance!

## Issues

Stability for Reynolds number larger than about 40.

More complex geometries.

I don't understand what is wall rotation, and what is pressure difference calculation.

## Acknowledgments

* Paper by the original author: http://dx.doi.org/10.13140/2.1.1606.1120