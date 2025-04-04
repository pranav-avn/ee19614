% Design of Proportional Integral (PI) Controller Parameters Using Ant Colony 
% Optimization Algorithm (ACO) for VMC Based Buck Converter 

% AIM: Designing a Proportional-Integral (PI) controller for a Voltage Mode Control (VMC) 
% based buck converter using Ant Colony Optimization (ACO) in MATLAB.

% Define buck converter parameters 
Vin = 10; % Input voltage 
Vout = 5; % Output voltage 
L = 1e-3; % Inductance 
C = 10e-6; % Capacitance 
R = 0.1; % Load resistance 

% Transfer function of the buck converter 
num = [Vout*L/R, Vout*(1-R*C/L), 0]; 
den = [L*C, R*C, 1]; 
G = tf(num, den); 

% Define ACO parameters 
num_ants = 10; % Number of ants 
num_iterations = 50; % Number of iterations 
evaporation_rate = 0.5; % Evaporation rate 
alpha = 1; % Pheromone importance factor 
beta = 1; % Heuristic information importance factor 

% Initialize pheromone matrix 
pheromones = ones(2, num_iterations); 
for iteration = 1:num_iterations 
    for ant = 1:num_ants 
        % Generate random values for PI controller parameters 
        Kp = rand; 
        Ki = rand; 
        % Simulate buck converter with current PI controller parameters 
        sys = feedback(G, Kp*tf([1, Ki], [1, 0])); 
        step_response = stepinfo(sys); 
             
        % Evaluate performance (e.g., overshoot, settling time) 
        performance = alpha * step_response.SettlingTime + beta * step_response.Overshoot; 
             
        % Update pheromone matrix 
        pheromones(:, iteration) = (1 - evaporation_rate) * pheromones(:, iteration) + performance; 
    end 
end 
 
% Extract best PI controller parameters 
[~, index] = min(pheromones); 
best_Kp = Kp(index); 
best_Ki = Ki(index); 
disp(['Best PI Controller Parameters: Kp = ', num2str(best_Kp), ', Ki = ', num2str(best_Ki)]); 

% Output:
% Best PI Controller Parameters: Kp = 0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019    0.060019
% 0.060019    0.060019    0.060019    0.060019    0.060019

% Ki = 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675     0.86675
% 0.86675     0.86675