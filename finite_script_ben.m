%% Exercise 4: Numerical Solution Using the Finite Difference Approach
% The following code implements the Finite Difference Method of solving
% Ordinary Differential Equations (ODE) using MATLAB(R)

%% Set Up
% Define ODE parameters, such as coefficients, f(t), time interval, etc.
%
% ODE Settings

clear all, close all        %Clears and closes all variables/figures

h = 0.01;                   %Width
A = 1;                      %Coefficient of x"
B = -20;                    %Coefficient of x'
C = 0;                      %Coefficient of x
t0 = 0; tf = 1;             %Interval
x0 = 1; xf = 0;             %Boundary Conditions 
func=@(t) -1;               %f(t) function handle
%% Calculate variables
% Constants determined from ODE Settings

t = [t0:h:tf];              %Time vector, increasing from t0 to tf by width h 
N = (tf-t0)/h;              %No of Iterations
%% 
% Evaluate _a_, _b_ and _c_, where
%
% $$ a = \frac{A}{h^2} - \frac{B}{2h} $$
% \indent $$ b = C - \frac{2A}{h^2} $$
% \indent $$ c = \frac{A}{h^2} + \frac{B}{2h} $$

a = (A/(h*h))-(B/(2*h));
b = C-((2*A)/(h*h));
c = (A/(h*h))+(B/(2*h));
%% 
% Assemble vector _vec_ with a *for loop*

vec = zeros(1,N-1);         %Vector initialised as array of 0s, size N-1
for i = 1:N-1
    vec(i) = func(t);       %Vector contents initialised to func(t) value
end

%% 
% Set first and last values of _vec_

vec(1) = func(t)-a*x0;
vec(N-1) = func(t)-c*xf;
%%
% Use *solvetridiag* function and add _x0_ and _xf_ to _y_ as _x_

y = solvetridiag(N-1,a,b,c,vec);
x = [x0;y;xf];
%% Create graph

figure(1)
plot(t,x, 'b*')
title('Solution using Finite Difference Method')
xlabel('Time t (s)'), ylabel('Output x')
