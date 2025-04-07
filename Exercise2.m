% Exercise 2
% Parameters
h = 0.125; 
t2 = 0:h:4; 
y_exact = cos(2*t2);

% Euler's Method 
y_euler = zeros(size(t2));
v_euler = zeros(size(t2)); 
y_euler(1) = 1; 
v_euler(1) = 0; 
for i = 2:length(t2)
    v_euler(i) = v_euler(i-1) - 4 * y_euler(i-1) * h;
    y_euler(i) = y_euler(i-1) + v_euler(i-1) * h;
end

% RK4 Method
y_rk4 = zeros(size(t2));
v_rk4 = zeros(size(t2));
y_rk4(1) = 1; 
v_rk4(1) = 0; 
for i = 2:length(t2)
    k1y = v_rk4(i-1);
    k1v = -4 * y_rk4(i-1);
    
    k2y = v_rk4(i-1) + 0.5 * h * k1v;
    k2v = -4 * (y_rk4(i-1) + 0.5 * h * k1y);
    
    k3y = v_rk4(i-1) + 0.5 * h * k2v;
    k3v = -4 * (y_rk4(i-1) + 0.5 * h * k2y);
    
    k4y = v_rk4(i-1) + h * k3v;
    k4v = -4 * (y_rk4(i-1) + h * k3y);
    
    v_rk4(i) = v_rk4(i-1) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);
    y_rk4(i) = y_rk4(i-1) + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
end

% Plotting
figure;
hold on;
plot(t2, y_euler, 'r-', 'DisplayName', "Euler's Method");
plot(t2, y_rk4, 'b-', 'DisplayName', "RK4 Method");
plot(t2, y_exact, 'k--*','DisplayName', "Exact Solution");
xlabel('Time (t)');
ylabel('y(t)');
legend;
title('Solutions of the Differential Equation');
grid on;
