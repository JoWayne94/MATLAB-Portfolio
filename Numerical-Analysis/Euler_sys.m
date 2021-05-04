clear; clc; close all

% Script to implement Explicit Euler to solve the model problem,  
% namely the initial value problem (IVP) of y' = lambda*y, y(0) = y_0 
% Written by Jo Wayne Tan, 20/01/2019. Academic responsible: Dr Thulasi Mylvaganam

width = 1.5 ; % Line thickness
lsize = 12 ; tsize = 12 ; legsize = 9 ; msize = 8 ; % Font size in plots, label, title, legend
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); set(0,'defaultLegendInterpreter','latex');

%% Problem definition
rho = 28;
sigma = 10;
beta = 8/3;
A = [-sigma sigma 0 ; rho -1 0 ; 0 0 -(beta)];
f = @(y) [sigma*(y(2)-y(1)); (y(1)*(rho-y(3))-y(2)); (y(1)*y(2))-((beta)*y(3))];
y0 = [1; 1; 1];  % The initial value y(0) 
tf = 50;         % The final time (up to which we want to integrate the ODE)

%% Parameters for the numerical solution
hmax = 2/max(abs(eig(A)));
h = 0.00001; % uniform time step size (delta t)
% h = hmax; 
% h = hmax + 0.01;

%% Numerical Solution (Explicit Euler) 
% The array y is the numerical solution of the equation. Initialise 
y_num1(:,1) = y0; % Initial condition. Note MATLAB starts indexing with `1` 
t_num1(1) = 0;    % Initial time (t = 0 in this case)
i = 1;            % Initialise counter

while t_num1(i) < tf - h 
     ydot = f(y_num1(:, i));                 % The is the ODE (y' = lambda*y) 
     y_num1(:, i+1) = y_num1(:, i) + h*ydot; % Euler method 
     t_num1(i+1) = t_num1(i) + h;            % Save "sample times" in an array
     i = i+1;   
end 

%% Plot numerical solution in red 
figure (1)
hold on; box on ; grid on 
plot3(t_num1, y_num1(1,:), y_num1(2,:), 'ro--', 'MarkerSize', msize, 'LineWidth', width)
plot(t_num1, y_num1(1,:), 'r', 'MarkerSize', msize, 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$y$', 'fontsize', lsize, 'Rotation', 0)
zlabel('$y_2$', 'fontsize', lsize, 'Rotation', 0)

subplot(2,1,1)
plot(t_num1, y_num1(2,:), 'r', 'MarkerSize', msize, 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$y_2$', 'fontsize', lsize, 'Rotation', 0)

subplot(2,1,2)
plot(t_num1, y_num1(3,:), 'r', 'MarkerSize', msize, 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$y_3$', 'fontsize', lsize, 'Rotation', 0)

%% Exact solution (known for the model problem) 
y1_sol = @(t) y0(1)*exp(A(1,1)*t);
y2_sol = @(t) y0(2)*exp(A(2,2)*t);
y3_sol = @(t) y0(3)*exp(A(3,3)*t);
N = 50000; 
t_exact = linspace(0, tf, N); % creates a vector of N time instants (uniformly spaced)
 
% Samples of the exact solution
y1_exact = y1_sol(t_exact); 
y2_exact = y2_sol(t_exact);
y3_exact = y3_sol(t_exact);
 
%% Plot exact solution (in blue) and compare with the numerical solution
figure (2)
plot3(y1_exact, y2_exact, y3_exact, 'b', 'LineWidth', width)

figure (3) 
hold on; box on; grid on
% subplot(2,1,1)
plot(t_exact, y1_exact, 'r', 'LineWidth', width)
% subplot(2,1,2)
plot(t_exact, y2_exact, 'b', 'LineWidth', width)
% subplot(2,1,3)
plot(t_exact, y3_exact, 'g', 'LineWidth', width)

%% 4th order Runge-Kutta
f = @(t,y) [sigma*(y(2)-y(1)); (y(1)*(rho-y(3))-y(2)); (y(1)*y(2))-((beta)*y(3))];
tspan = [0, tf];
y0 = [1; 1; 1];
h = 0.01; % uniform time step size (delta t)

m = length(y0);
t = tspan(1):h:tspan(2);
y = zeros(m, length(t));
y(:,1) = y0;

for i = 1:length(t) - 1
    k1 = h.*f(t(i), y(:,i));
    k2 = h.*f(t(i)+h/2, y(:,i) + (k1/2));
    k3 = h.*f(t(i)+h/2, y(:,i) + (k2/2));
    k4 = h.*f(t(i)+h, y(:,i) + k3);
    y(:,i+1) = y(:,i) + k1/6 + (k2 + k3)/3 + k4/6;
end

figure (4)
subplot(3,1,1)
hold on; box on; grid on
plot(t_num1, y_num1(1,:), 'r', 'LineWidth', width)
plot(t, y(1,:), 'b', 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$x$', 'fontsize', lsize, 'Rotation', 0)

subplot(3,1,2)
hold on; box on; grid on
plot(t_num1, y_num1(2,:), 'r', 'LineWidth', width)
plot(t, y(2,:), 'b', 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$y$', 'fontsize', lsize, 'Rotation', 0)

subplot(3,1,3)
hold on; box on; grid on
plot(t_num1, y_num1(3,:), 'r', 'LineWidth', width)
plot(t, y(3,:), 'b', 'LineWidth', width)
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('$z$', 'fontsize', lsize, 'Rotation', 0)

figure (5)
hold on; box off; grid on 
xlabel('x(t)', 'fontsize', lsize)
ylabel('y(t)', 'fontsize', lsize)
zlabel('z(t)', 'fontsize', lsize)
plot3(y(1,:), y(2,:), y(3,:), 'b', 'LineWidth', width) 
plot3(y_num1(1,:), y_num1(2,:), y_num1(3,:), 'r', 'LineWidth', width)
