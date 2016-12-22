# Particle Swarm Optimization

This repository contains the implementation of several Particle Swarm Optimization (PSO) techniques in Octave.

## Implementation
PSO helps to find solutions for a wide range of problems and works without traditional optimization methods such as gradient descent/ascent.
In order to show the effectiveness of PSO, we decided to create two different scenarios. The first is a classic minimization/maximization problem, while the second implementation will try to solve the travelling salesman problem (TSP) using a fuzzy probabilistic algorithm.

### Simple maximization/minimization (PSO1)
The file PSO1.m contains a basic implementation of the core PSO algorithm. To test it, we decided to minimize two different functions:
* the [Sphere function](https://www.sfu.ca/~ssurjano/spheref.html)
* the [Rastrigin Function](https://www.sfu.ca/~ssurjano/rastr.html)

The parameters used can be seen in the Octave file. Our implementation manages to minimize both functions rather fast. In order to test the algorithm with more dimensions, the *nr_variables* variable can be adjusted.

### More complex PSO (Travelling Salesman Problem)
For our more elevated implementation we decided to solve (or at least approximate) the travelling salesman problem (TSP). Our implementation is based on a modification of the original PSO algorithm, which uses a fuzzy position matrix and probabilities in order to choose a route. More information about the core idea can be found in [1] and [2].
Our implementation is not always able to find the optimal solution, but it approximates it quite well, given that the PSO algorithm was originally not build to solve TSP, an NP-hard problem.
We used different datasets to test our implementation, mainly three standard ones:
- [burma14](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/burma14.tsp)
- [ulysses](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/ulysses.tsp)
- [berlin52](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/berlin52.tsp)

For ulysses there are more then a trillion possible routes, and yet PSO manages to approximate it very well in just a few seconds.

## PSO Termination Conditions
There are several conditions that can be used for the PSO algorithm in order to know when to stop:
- Maximum iteration number: The PSO algorithm stops after a set number of iterations.
- Number of iterations without change: The PSO algorithm stops after a set amount of iterations, which didn't change the result.
- Minimum cost function difference: The PSO algorithm stops after the difference between the last obtained cost and the current best cost is less then a set threshold.

We have chosen to use a fixed maximum number of iterations for our implementations.
