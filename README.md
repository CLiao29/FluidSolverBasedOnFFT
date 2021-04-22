# FluidSolverBasedOnFFT
This project is based on

[*J. Stam, A Simple Fluid Solver based on the FFT, Journal of Graphics Tools
Volume 6, Number 2, 2001, 43-52*](www.baidu.com)


**Abstract** 

The goal of this project is to introduce a fluid solver based on the Fast Fourier Transform (FFT). The main algorithm of the solver, similar to many other fluid solvers, are built based on the Naive-Stokes equations, which are a set of partial differential equations describing the motion of viscous fluid simulations. The fluids are represented in a n*n computational grids, with each grid having the corresponding horizontal and vertical velocity vectors. The main advantage is that, all the computational process can be made simple and concise, and are concatenated together into one method. Any user who want to make the simulation can iterate through the solver method again and again with the initial external force provided. The project is originally implemented in C and it utilizes the the MIT's FFTW for the solver of the fast fourier transform. In this project that I presented, the FFT solver is based on the Jtransform external jar in Java. Besides, a Complex class is used to represent the abstraction of the complex number. A simple GUI, showing the actual vector fields of each grid, is also included. 
