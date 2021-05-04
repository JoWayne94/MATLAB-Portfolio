clear; clc; close all

% Script to implement Explicit Euler to solve the model problem,  
% namely the initial value problem (IVP) of y' = lambda*y, y(0) = y_0 
% Written by Jo Wayne Tan, 20/01/2019. Academic responsible: Dr Thulasi Mylvaganam

width = 1.5 ; % Line thickness
lsize = 12 ; tsize = 12 ; legsize = 9 ; msize = 8 ; % Font size in plots, label, title, legend
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); set(0,'defaultLegendInterpreter','latex');

%% Problem definition
lambda = -1;
y0 = 1;  % The initial value y(0) 
tf = 8;  % The final time (up to which we want to integrate the ODE)

%% Parameters for the numerical solution
% The parameters (time step, etc) for the numerical solution can be
% varied to explore the effects on accuracy and numerical stability  

h = 0.1; % uniform time step size, delta t
% To compare different time steps simultaneously, you may wish to add
% alternative time steps h2, h3, etc here and repeat the numerical solution
% loop for these (in the next block). Note you will have to give t and y a
% a different name to do this.

%% Numerical Solution (Explicir Euler) 
% The array y is the numerical solution of the equation
% Initialise 
y_num1(1) = y0; % Initial condition. Note MATLAB starts indexing with one 
t_num1(1) = 0;  % Initial time (t = 0 in this case) 
i = 1;          % Initialise counter 

while t_num1(i) <= tf 
     ydot = lambda * y_num1(i);        % The is the ODE (y' = lambda*y) 
     y_num1(i+1) = y_num1(i) + h*ydot; % Euler method 
     t_num1(i+1) = t_num1(i) + h;      % Save "sample times" in an array 
     
     i = i+1;     
end 

% Plot numerical solution in red 
figure(1) 
hold on ; box on ; grid on
plot1 = plot(t_num1, y_num1, 'ro--','MarkerSize', msize, 'LineWidth', width);
xlabel('time (s)', 'fontsize', lsize, 'Rotation', 0)
ylabel('y', 'fontsize', lsize, 'Rotation', 0)
title('Numerical vs Exact Solution of $\dot{y} = \lambda * y$ using Explicit Euler', 'FontSize', tsize);

%% Exact Solution
% The exact solution (known for the model problem) 
y_sol = @(t) y0 * exp(lambda * t); % exact solution of the problem 
N = 500;
t_exact = linspace(0, tf, N); % creates a vector of N time instants (uniformly spaced) 
y_exact = y_sol(t_exact);     % Samples of the exact solution  
 
% Plot exact solution (in blue) to compare with the numerical solution
figure(1) 
plot2 = plot(t_exact, y_exact, 'b', 'LineWidth', width);
L = legend([plot1 plot2], 'Numerical Soln', 'Exact Soln', 'Location', 'best');
L.FontSize = legsize;