# Particle Swarm Optimization

This repository contains the implementation of several Particle Swarm Optimization (PSO) techniques in Octave.

## Implementation
PSO helps to find solutions for a wide range of problems and works without traditional optimization methods such as gradient descent/ascent.
In order to show the effectiveness of PSO, we decided to create two different scenarios. The first is a classic minimization/maximization problem, while the second scenario is a more complex one, specifically tailored for PSO.

## PSO Termination Conditions
There are several conditions that can be used for the PSO algorithm in order to know when to stop:
- Maximum iteration number: The PSO algorithm stops after a set number of iterations.
- Number of iterations without change: The PSO algorithm stops after a set amount of iterations, which didn't change the result.
- Minimum cost function difference: The PSO algorithm stops after the difference between the last obtained cost and the current best cost is less then a set threshold.

We have chosen to use a fixed maximum number of iterations for our implementations.
