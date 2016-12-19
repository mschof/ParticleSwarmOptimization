function PSO1
  
  %% Clear all
  clc;
  clear;
  close all;
  
  %% Problem Definition
  cf = @RastriginFunction;
  nr_variables = 5;                       % Number of variables unknown (part of the decision)
  variable_size = [1 nr_variables];        % Vector representation
  var_min = -10;                          % Lower bound of decision space
  var_max = 10;                           % Upper bound of decision space
  
  %% Parameter Adjustment
  max_iterations = 100;                   % Maximum iterations in PSO algorithm
  swarm_size = 50;                        % Swarm size (number of particles)
  w = 1;                                  % Inertia coefficient
  w_damp = 0.90;                          % damping of inertia coefficient, lower = faster damping
  c1 = 1;                                 % Cognitive acceleration coefficient (c1 + c2 = 4)
  c2 = 3;                                 % Social acceleration coefficient (c1 + c2 = 4)
  
  %% Init
  template_particle.position = [];
  template_particle.velocity = [];
  template_particle.cost = [];
  template_particle.best.position = [];   % Local best
  template_particle.best.cost = inf;       % Local best
  
  % Copy and put the particle into a matrix
  particles = repmat(template_particle, swarm_size, 1);
  
  % Initialize global best (current worst value, inf for minimization, -inf for maximization)
  global_best.cost = inf;
  
  for i=1:swarm_size
    
    % Initialize all particles with random position inside the search space
    particles(i).position = unifrnd(var_min, var_max, variable_size);
    
    % Initiliaze velocity to the 0 vector
    particles(i).velocity = zeros(variable_size);
    
    % Evaluate the current cost
    particles(i).cost = cf(particles(i).position);
    
    % Update the local best to the current location
    particles(i).best.position = particles(i).position;
    particles(i).best.cost = particles(i).cost;
    
    % Update global best
    if (particles(i).best.cost < global_best.cost)
      global_best.position = particles(i).best.position;
      global_best.cost = particles(i).best.cost;
    endif
    
  endfor
  
  % Best cost at each iteration
  best_costs = zeros(max_iterations, 1);
    
  
  %% PSO Loop
  
  for iteration=1:max_iterations
    
    for i=1:swarm_size
      
      % Initialize two random vectors
      r1 = rand(variable_size);
      r2 = rand(variable_size);
      
      % Update velocity
      particles(i).velocity = (w * particles(i).velocity) ...
        + (c1 * r1 .* (particles(i).best.position - particles(i).position)) ...
        + (c2 * r2 .* (global_best.position - particles(i).position));
      
      % Update position
      particles(i).position = particles(i).position + particles(i).velocity;
      
      % Update cost
      particles(i).cost = cf(particles(i).position);
      
      % Update local best (and maybe global best) if current cost is better
      if (particles(i).cost < particles(i).best.cost)
        particles(i).best.position = particles(i).position;
        particles(i).best.cost = particles(i).cost;
        
        % Update global best
        if (particles(i).best.cost < global_best.cost)
          global_best.position = particles(i).best.position;
          global_best.cost = particles(i).best.cost;
        endif
        
      endif
      
    endfor
    
    % Get best value
    best_costs(iteration) = global_best.cost;
    
    % Display information for this iteration
    % disp(["Iteration " num2str(iteration) ": best cost = " num2str(best_costs(iteration))]);
    
    % Damp w
    w = w * w_damp;
    
  endfor
  
  
  %% Print Results
  figure;
  plot(best_costs, "LineWidth", 2);
  xlabel("iteration");
  ylabel("best cost");
  
endfunction

function ret = SphereFunction(x)
  
  ret = sum(x.^2);

endfunction

function ret = RastriginFunction(x)
  
  ret = 10 * size(x)(2) + sum(x.^2 - 10 * cos(2*pi.*x));

endfunction
