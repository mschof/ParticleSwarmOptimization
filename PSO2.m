% Remarks: Some TSP problems should use the geographic distance rather than the euclidean one!

function PSO2

  %% Clear all
  clc;
  clear;
  close all;

  % Import TSP data
  tsp_data;

  %% Problem Definition
  % Coordinates of cities (see also: http://www.cc.gatech.edu/~bdilkina/CSE6140-2014fa/ProjectTSP.htm)
  data = berlin52;%burma14;%berlin52;%%ulysses16;
  coords = data.coords;
  mode = data.mode;
  nr_coords = length(coords);

%  testprobs = ...
%    [0 1 0 0 0 0 0 0 0 0 0 0 0 0;
%    0 0 1 0 0 0 0 0 0 0 0 0 0 0;
%    0 0 0 1 0 0 0 0 0 0 0 0 0 0;
%    0 0 0 0 1 0 0 0 0 0 0 0 0 0;
%    0 0 0 0 0 1 0 0 0 0 0 0 0 0;
%    0 0 0 0 0 0 1 0 0 0 0 0 0 0;
%    0 0 0 0 0 0 0 1 0 0 0 0 0 0;
%    0 0 0 0 0 0 0 0 1 0 0 0 0 0;
%    0 0 0 0 0 0 0 0 0 1 0 0 0 0;
%    0 0 0 0 0 0 0 0 0 0 1 0 0 0;
%    0 0 0 0 0 0 0 0 0 0 0 1 0 0;
%    0 0 0 0 0 0 0 0 0 0 0 0 1 0;
%    0 0 0 0 0 0 0 0 0 0 0 0 0 1;
%    1 0 0 0 0 0 0 0 0 0 0 0 0 0;]
%  testroute = calcroute(testprobs)
%  cost = TSPCost(coords, testroute)

%  testprobs = initposition(nr_coords)
%  testroute = calcroute(testprobs)
%  cost = TSPCost(coords, testroute)

  cf = @TSPCost;

  %% Parameter Adjustment (adjust them for better performance!)
  max_iterations = 500;                     % Maximum iterations in PSO algorithm (use 1000 later)
  swarm_size = 50;                        % Swarm size (number of particles)
  w_high = 1.0;
  w_low = 0.2;
  w = w_low;                             % Inertia coefficient, TODO: find the best value here!
  w_damp = 1.0;                          % damping of inertia coefficient, lower = faster damping, TODO: find the best value here!
  c1 = 2;                                 % Cognitive acceleration coefficient (here: 0 >= c1 <= 2)
  c2 = 2;                                 % Social acceleration coefficient (here: 0 >= c2 <= 2)

  % Simulated Annealing Params
  T_init = 1000;
  T_end = 1;
  T_rate = 0.99;
  T_cur = T_init;
  SA_m = 0;
  SA_gen = 30;

  %% Init
  template_particle.position = zeros(nr_coords);
  template_particle.velocity = zeros(nr_coords);
  template_particle.cost = 0;
  template_particle.best.position = zeros(nr_coords);   % Local best
  template_particle.best.cost = inf;                    % Local best

  % Copy and put the particles into a matrix
  particles = repmat(template_particle, swarm_size, 1);

  % Initialize global best (current worst value, inf for minimization, -inf for maximization)
  global_best.cost = inf;
  ultimate_best.cost = inf;

  % random search
  % cost_t = cf(coords, calcroute(initposition(nr_coords)), mode);
  % while (cost_t > 16000)
  %   cost_t = cf(coords, calcroute(initposition(nr_coords)), mode);
  % endwhile
  % return;

  for i=1:swarm_size

    % Initialize all particles with random position inside the search space
    particles(i).position = initposition(nr_coords);

    % Initiliaze velocity to the 0 vector
    particles(i).velocity = initvelocity(nr_coords);

    % Evaluate the current cost
    particles(i).cost = cf(coords, calcroute(particles(i).position), mode);

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
    global_best_new.position = global_best.position;
    global_best_new.cost = global_best.cost;
    for i=1:swarm_size

      % Initialize two random vectors
      r1 = rand();
      r2 = rand();

      % Update velocity
      particles(i).velocity = (w * particles(i).velocity) ...
        + (c1 * r1 .* (particles(i).best.position - particles(i).position)) ...
        + (c2 * r2 .* (global_best.position - particles(i).position));

%      % should all be 0...
%      round(sum((w * particles(i).velocity)(1, :)))
%      round(sum((particles(i).best.position - particles(i).position)(1, :)))
%      round(sum((global_best.position - particles(i).position)(1, :)))
%      round(sum(particles(i).velocity(1, :)))

      % Update position
      particles(i).position = particles(i).position + particles(i).velocity;

      % Normalize position
      particles(i).position = normprobabilities(particles(i).position);

      % Update cost
      particles(i).cost = cf(coords, calcroute(particles(i).position), mode);

      % Update local best (and maybe global best) if current cost is better
      if (particles(i).cost < particles(i).best.cost)
        particles(i).best.position = particles(i).position;
        particles(i).best.cost = particles(i).cost;

        % Update new global best
        if (particles(i).best.cost < global_best.cost)
          global_best_new.position = particles(i).best.position;
          global_best_new.cost = particles(i).best.cost;
        endif

      endif

    endfor

    % Update global best
    if (global_best_new.cost < global_best.cost)
      global_best.position = global_best_new.position;
      global_best.cost = global_best_new.cost;
    else
      % No improvement...
      SA_m = SA_m + 1;
    endif

    % Update ultimate best
    if (global_best.cost < ultimate_best.cost)
      ultimate_best.position = global_best.position;
      ultimate_best.cost = global_best.cost;
    endif

    % Get best value
    best_costs(iteration) = global_best.cost;

    % Simulated Annealing
    if (SA_m == SA_gen)
      % Annealing...
      while (T_cur > T_end)
        neighbour_solution = createneighbour(global_best.position);
        neighbour_cost = cf(coords, calcroute(neighbour_solution), mode);
        difference = neighbour_cost - global_best.cost;
        if (min(1, exp(-difference / T_cur)) > rand())
          global_best.position = neighbour_solution;
          global_best.cost = neighbour_cost;
        endif
        T_cur = T_cur * T_rate;
      endwhile
      SA_m = 0;
    endif

    % Display information for this iteration
    % disp(["Iteration " num2str(iteration) ": best cost = " num2str(best_costs(iteration))]);

    % Damp w
    w = w * w_damp;
    %w = w - ((w_high - w_low) / max_iterations);

  endfor

  %% Print results
  bestroute = calcroute(ultimate_best.position);
  ["Best cost: " num2str(ultimate_best.cost)]
  ["Route: " num2str(bestroute)]
  ["Error: " num2str((((ultimate_best.cost / data.optimum) - 1) * 100)) "%"]
  plotroute(coords, bestroute);

  %% Plot results
  figure;
  plot(best_costs, "LineWidth", 2);
  xlabel("iteration");
  ylabel("best cost");

endfunction

function ret = initposition(nr_coords)

  pos = zeros(nr_coords);

  for i=1:nr_coords

    row_rand = rand(1, nr_coords - 1); % diagonal should be very small (-inf), so a city can't lead to itself in the next step, therefore (nr_coords - 1) random values
    row_rand = row_rand / sum(row_rand); % normalize so that sum = 1
    pos(i, 1:(i - 1)) = row_rand(1:(i - 1));
    pos(i, i) = 0; % TODO: 0 instead of -inf?
    pos(i, (i + 1):nr_coords) = row_rand(i:(nr_coords - 1));

  endfor

  ret = pos;

endfunction

function ret = initvelocity(nr_coords)

  pos = zeros(nr_coords); % is this alone maybe also enough?

  for i=1:nr_coords

    for j=1:(nr_coords / 2)

      r1 = rand();
      r2 = r1 * -1.0;
      pos(i, j) = r1;
      pos(i, nr_coords - (j - 1)) = r2;

    endfor

    % Make center value 0 if nr_coords is odd
    if (mod(nr_coords, 2) == 0)
      pos(i, (nr_coords / 2) + 1) = 0;
    endif

    % Shuffle
    pos(i, :) = pos(i, :)(randperm(nr_coords));

  endfor

  ret = pos;

endfunction

function ret = createneighbour(probabilities)

  neighbour_solution = probabilities;
  r_index_1 = randi([1, length(probabilities)]);
  r_index_2 = randi([1, length(probabilities)]);
  while (r_index_1 == r_index_2)
    r_index_2 = randi([1, length(probabilities)]);
  endwhile
  neighbour_solution([r_index_1, r_index_2], :) = neighbour_solution([r_index_2, r_index_1], :);
  ret = neighbour_solution;

endfunction

function ret = calcroute(probabilities)

  % I think that one should take the max of all columns *which are not already taken* in order to create a valid route every time
  nr_rows = size(probabilities, 1); % This should correspond to the number of coords
  route = zeros(1, nr_rows + 1);
  route(1) = 1; % The first city is always the first city :D

  for i=1:(nr_rows - 1) % TODO: Is it also important to regard the last row? Because the final route (end -> begin) should be clear...

    current_row = probabilities(i, :);
    current_row(route(route ~= 0)) = -1; % set used indices to -1
    [max_element, max_index] = max(current_row); % max_index corresponds to the index of the city to "move next"
    route(i + 1) = max_index;

  endfor

  % Add begin to route again (end -> begin as last step)
  route(length(route)) = 1;

  ret = route;

endfunction

function ret = normprobabilities(probabilities)

  % convert negative numbers to 0
  probabilities(probabilities < 0) = 0;

  % normalize resulting rows so that the sum of each row = 1
  for i=1:length(probabilities)

    s = sum(probabilities(i, :));
    if (s != 0)
      probabilities(i, :) = probabilities(i, :) / sum(probabilities(i, :));
    endif

  endfor

  ret = probabilities;

endfunction

function ret = geodist(p1, p2)
  % function from http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/TSPFAQ.html

  deg_p1 = floor(p1);
  minute_p1 = p1 - deg_p1;
  rad_p1 = pi * (deg_p1 + 5.0 * minute_p1 / 3.0) / 180.0;
  deg_p2 = floor(p2);
  minute_p2 = p2 - deg_p2;
  rad_p2 = pi * (deg_p2 + 5.0 * minute_p2 / 3.0) / 180.0;

  earth_radius = 6378.388; % Hayford ellipsoid

  q1 = cos(rad_p1(2) - rad_p2(2));
  q2 = cos(rad_p1(1) - rad_p2(1));
  q3 = cos(rad_p1(1) + rad_p2(1));
  distance = floor((earth_radius * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0));

  ret = distance;

endfunction

function ret = TSPCost(coords, route, mode)

  total_cost = 0;
  % Add distances of route (begin -> end is already included, so first entry and last entry should both be 1)
  for i=1:(length(route) - 1)
    current = [coords(route(i), :)(2) coords(route(i), :)(3)];
    next = [coords(route(i + 1), :)(2) coords(route(i + 1), :)(3)];
    distance_current = 0;
    if (strcmp(mode, "eucl")) % euclidean distance
      distance_current = norm(next - current);
    elseif (strcmp(mode, "geo")) % geographical distance using latitude/longitude and Hayford ellipsoid
      distance_current = geodist(next, current);
    else
      ["Invalid distance calculation mode " mode]
      ret = -1;
      return;
    endif
    total_cost = total_cost + distance_current;
  endfor

  ret = total_cost;

endfunction

function plotroute(coords, route)

  edgesX = zeros(1, length(route));
  edgesY = zeros(1, length(route));
  for i=1:length(route)

    edgesX(i) = coords(route(i), 2);
    edgesY(i) = coords(route(i), 3);

  endfor

  plot(edgesX, edgesY);

endfunction
