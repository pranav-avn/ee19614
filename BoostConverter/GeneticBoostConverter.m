% Design of Proportional Integral (PI) Controller Parameters Using Genetic 
% Algorithm (GA) for VMC Based Boost Converter 

% AIM: Designing a Proportional Integral (PI) controller for a Variable Mode Control (VMC) 
% based boost converter using a Genetic Algorithm (GA).

% Requires: Global Optimisation Toolbox

% Define boost converter parameters 
Vin = 12;   % Input voltage 
Vout = 24;  % Output voltage       
L = 100e-6; % Inductance (H)     
C = 100e-6; % Capacitance (F)     
R = 0.1;    % Load resistance (ohm)       
fs = 50e3;  % Switching frequency (Hz)       
D = 0.5;    % Initial duty cycle        
T = 1/fs;   % Switching period    

% Define simulation time 
t = 0:T:1/fs*10; % Simulate for 10 switching cycles 

% Objective function to minimize the steady-state error 
objectiveFunction = @(K) boostConverterCostFunction(Vin, Vout, L, C, R, fs, t, K(1), K(2)); 

% Define GA parameters 
populationSize = 50; 
numberOfGenerations = 100; 
mutationRate = 0.1; 

% Lower and upper bounds for Kp and Ki 
lb = [0 0]; 
ub = [10 10]; 

% Run GA to optimize PI controller parameters 
options = optimoptions('ga', 'PopulationSize', populationSize, 'MaxGenerations', numberOfGenerations, 'MutationFcn', {@mutationadaptfeasible, mutationRate}); 
[K_opt, cost_opt, exitflag, output] = ga(objectiveFunction, 2, [], [], [], [], lb, ub, [], options); 

% Display optimization results 
if exitflag > 0 
    disp('Optimization successful.'); 
    disp('Optimized PI controller parameters (Kp, Ki):'); 
    disp(K_opt); 
    disp('Optimized cost:'); 
    disp(cost_opt); 
else 
    disp('Optimization failed to converge.'); 
end 

% Simulate boost converter with optimized PI controller parameters 
[vout, il, vc] = boostConverterSimulation(Vin, Vout, L, C, R, fs, t, D, K_opt(1), K_opt(2)); 

% Plot results 
figure; 
subplot(3,1,1); 
plot(t, vout); 
xlabel('Time (s)'); 
ylabel('Output Voltage (V)'); 
title('Boost Converter Output Voltage'); 

subplot(3,1,2); 
plot(t, il); 
xlabel('Time (s)'); 
ylabel('Inductor Current (A)'); 
title('Boost Converter Inductor Current'); 

subplot(3,1,3); 
plot(t, vc); 
xlabel('Time (s)'); 
ylabel('Capacitor Voltage (V)'); 
title('Boost Converter Capacitor Voltage'); 

% Function to calculate steady-state error as cost function 
function cost = boostConverterCostFunction(Vin, Vout, L, C, R, fs, t, Kp, Ki) 
    [~, ~, vc] = boostConverterSimulation(Vin, Vout, L, C, R, fs, t, 0.5, Kp, Ki); % Simulate with a fixed initial duty cycle 
    error = Vout - vc(end); % Steady-state error 
    cost = error^2; % Squared error as cost function 
end 

% Function to simulate boost converter operation 
function [vout, il, vc] = boostConverterSimulation(Vin, Vout, L, C, R, fs, t, D, Kp, Ki) 
    vout = zeros(size(t)); 
    il = zeros(size(t)); 
    vc = zeros(size(t)); 
    for i = 1:length(t) 
        if mod(t(i), 1/fs) < D*(1/fs) 
            duty_cycle = 1; 
        else 
            duty_cycle = 0; 
        end 
        if i == 1 
            il_prev = 0; 
            vc_prev = 0; 
        else 
            il_prev = il(i-1); 
            vc_prev = vc(i-1); 
        end 
        if duty_cycle == 1 
            il(i) = il_prev + ((Vin - Vout) / L) * (1/fs); 
        else 
            il(i) = il_prev - (Vout / L) * (1/fs); 
    end 
    if duty_cycle == 1 
        vc(i) = vc_prev - (Vout / C) * (1 - exp(-1/(R * C) * (1/fs))) * (1/fs); 
    else 
        vc(i) = vc_prev + ((Vin - Vout) / C) * (1 - exp(-1/(R * C) * (1/fs))) * (1/fs); 
    end 
    vout(i) = duty_cycle * Vin - (1 - duty_cycle) * vc(i) - il(i) * R; 
    end 
end

% Output
% GeneticBoostConverter
% Optimization successful.
% Optimized PI controller parameters (Kp, Ki):
%          0    0.8404
%
% Optimized cost:
%    4.8517e+03