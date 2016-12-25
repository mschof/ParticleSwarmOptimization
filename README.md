# Particle Swarm Optimization

This repository contains the implementation of several Particle Swarm Optimization (PSO) techniques in Octave.

## Implementation
PSO helps to find solutions for a wide range of problems and works without traditional optimization methods such as gradient descent/ascent.
In order to show the effectiveness of PSO, we decided to create two different scenarios. The first is a classic minimization/maximization problem, while the second implementation will try to solve the travelling salesman problem (TSP) using a fuzzy probabilistic algorithm.

### Simple maximization/minimization (PSO1)
The file PSO1.m contains a basic implementation of the core PSO algorithm. To test it, we decided to minimize two different functions:
* [Sphere function](https://www.sfu.ca/~ssurjano/spheref.html)
* [Rastrigin Function](https://www.sfu.ca/~ssurjano/rastr.html)

The parameters used can be seen in the Octave file. Our implementation manages to minimize both functions rather fast. In order to test the algorithm with more dimensions, the *nr_variables* variable can be adjusted.

### More complex PSO with an approximation for the Travelling Salesman Problem (PSO2)
For our more elevated implementation we decided to solve (or at least approximate) the travelling salesman problem (TSP). Our implementation is based on a modification of the original PSO algorithm, which uses a fuzzy position matrix and probabilities in order to choose a route. More information about the core idea can be found in [1] and [2]. We have also implemented the [Simulated Annealing](https://en.wikipedia.org/wiki/Simulated_annealing) addition, because it leads to much better results.

Our implementation is not always able to find the optimal solution, but it approximates it quite well, given that the PSO algorithm was originally not build to solve TSP, an NP-hard problem.
We used different datasets to test our implementation, mainly three standard ones:
- [burma14](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/burma14.tsp)
- [ulysses16](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/ulysses16.tsp)
- [berlin52](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/berlin52.tsp)

We assumed the [Hayford ellipsoid](https://en.wikipedia.org/wiki/Hayford_ellipsoid) for all geographic problems.
For *ulysses16* there are more than a trillion possible routes, and yet PSO manages to approximate it quite good in a relatively low amount of time.

## PSO Termination Conditions
There are several conditions that can be used for the PSO algorithm in order to know when to stop:
- Maximum iteration number: The PSO algorithm stops after a set number of iterations.
- Number of iterations without change: The PSO algorithm stops after a set amount of iterations, which didn't change the result.
- Minimum cost function difference: The PSO algorithm stops after the difference between the last obtained cost and the current best cost is less then a set threshold.

We have chosen to use a fixed maximum number of iterations for our implementations.

## References
[1] W. Pang, K.-P. Wang, C.-G. Zhou and L.-J. Dong, "Fuzzy Discrete Particle Swarm Optimization for Solving Traveling Salesman Problem" ([Link](https://ai2-s2-pdfs.s3.amazonaws.com/3a30/480f7ccccecb02dbd951fe217eb64db5cd66.pdf))

[2] R. F. Abdel-Kader, "Fuzzy Particle Swarm Optimization with Simulated Annealing and Neighborhood Information Communication for Solving TSP" ([Link](http://thesai.org/Downloads/Volume2No5/Paper%203-Fuzzy%20Particle%20Swarm%20Optimization%20with%20Simulated%20Annealing%20and%20Neighborhood%20Information%20Communication%20for%20Solving%20TSP.pdf))
