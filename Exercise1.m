% Exercise 1
% Parameters
kgm = 0.026; 
pmax = 12000; 
p0 = 2555; 
t_data = 1950:10:2000; 
p_data = [2555, 3040, 3708, 4454, 5276, 6079]; 
step = 1; 
t = 1950:step:2000;

% Euler's Method
p_euler = zeros(size(t));
p_euler(1) = p0;
for i = 2:length(t)
    dpdt = kgm * (1 - p_euler(i-1) / pmax) * p_euler(i-1);
    p_euler(i) = p_euler(i-1) + dpdt * step;
end

% Heun's Method
p_heun = zeros(size(t));
p_heun(1) = p0;
for i = 2:length(t)
    dpdt1 = kgm * (1 - p_heun(i-1) / pmax) * p_heun(i-1);
    p_predictor = p_heun(i-1) + dpdt1 * step;
    dpdt2 = kgm * (1 - p_predictor / pmax) * p_predictor;
    p_heun(i) = p_heun(i-1) + (dpdt1 + dpdt2) / 2 * step;
end

% Plotting
figure;
hold on;
plot(t, p_euler, 'r-', 'DisplayName', "Euler's Method");
plot(t, p_heun, 'b-', 'DisplayName', "Heun's Method");
plot(t_data, p_data, 'ko-', 'DisplayName', "Measured Data");
xlabel('Year');
ylabel('Population/million');
legend;
title('Population Simulation');
grid on;
