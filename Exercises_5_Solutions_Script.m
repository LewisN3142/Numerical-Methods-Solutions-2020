% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long;

%% Class 5 Exercise 1 Script
% This script takes the Butcher Tableaux given in the question and uses the
% corresponding RK method to solve the ODEs (i) and (ii), plotting the
% first against its analytical solution, and the phase portrait for the
% second.

% Butcher Tableaux from question.
A = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0];
b = [1/6;1/3;1/3;1/6];
c = [0;1/2;1/2;1];

% Question (i)

% Exact solution
xExact = @(t) 1./(t.^2 +1); 
% RHS of ODE
f = @(t,x) -2*t*x.^2;

% RK with initial value (0,1), step size 0.5, and 100 steps
[t1,x1] = C5Ex1_RKfull(f,1,0,0.5,100,A,b,c);

figure(1);clf;hold on;
plot(t1,x1)
plot(t1,xExact(t1))
legend('RKapprox','Analytic')
title('Plot of Analytic and RK solutions of ODE(i)')
xlabel('t')
ylabel('x(t)')

figure(2);clf;hold on;
plot(xExact(t1),x1)
title('Plot of Analytic against RK solution of ODE(i)')
xlabel('Analytic solution')
ylabel('RK Approximation')

figure(3);clf;hold on;
plot(t1,xExact(t1)-x1)
title('Error in approximate solution')
xlabel('t')
ylabel('error')


% Question (ii)

% RHS of ODE (system of 2 equations so vector)
f = @(t, x) [0, 1; -1, 0]*x;

% RK with initial value (0, [1;0]), step size 0.5, and 1000 steps
[t2,x2] = C5Ex1_RKfull(f,[1;0],0,0.1,1000,A,b,c);

% Plot the components of x against t (displays cos and sine)
figure(4);hold on;
plot(t2,x2) 
legend('x component','y component')
title('Plot of components of RK solution of ODE(ii)')
xlabel('t')
ylabel('x(t) components')

% Plot the vector x in the cartesian plane as phase portrait (displays circle)
figure(5);hold on;
plot(x2(1,:),x2(2,:)) 
title('Plot of RK solution of ODE(ii) as phase portrait')
xlabel('x(1)')
ylabel('x(2)')


%% Class 5 Exercise 2 Script
% Maths exercise, see pdf for solution.


%% Class 5 Exercise 3 Script
% Execute this script to view stability plots for RK4 and Euler

C5Ex3_plotStability();


%% Class 5 Exercise 1 Functions
% A series of functions to perform Runge--Kutta for a given ODE and Butcher
% Tableau.

% This function carries out the full Runge--Kutta process for initial data (t0,x0), 
% Butcher tableau (A,b,c), function f, step size h, and number of steps m.%
% By increasing the number of steps we get an approximation on a larger interval
% By decreasing the step size, we get a higher resolution and more precise approx

function [t,x] = C5Ex1_RKfull(f,x0,t0,h,numSteps,A,b,c)
    
% Vector containing time-steps for the integration
t = t0 + h*(0:numSteps);
% Matrix where column i stores value of x after i-1 steps of RK
x = zeros(length(x0),numSteps+1);   
x(:,1) = x0;

% m steps of RK
for i = 1:numSteps
    x(:,i+1) = C5Ex1_RKstep(f,x(:,i),t(i),h,A,b,c); 
end

end 


% This function computes a single step of explicit RK for a function f,
% initial data (t0,x0), Butcher tableau (A,b,c) and step size h.
% f represents the RHS of a system of ODEs.

function x1 = C5Ex1_RKstep(f, x0, t0, h, A, b, c)

% First check if Butcher tableau was input correctly.
if ~iscolumn(b) 
    error('Sixth input b must be column vector')    
end
if ~iscolumn(c)
    error('Last input c must be column vector')    
end
n = size(A,[1,2]);
if n(1) ~= n(2)
    error('Fifth input A must be square matrix')
end

    % Number of variables/ equations
    n = length(x0);
    % Order of algorithm
    s = length(b);
    % Matrix, columns will be the k_i, will be d of these
    k = zeros(n, s);

    %Populate the matrix k
    for i = 1:s
    % Here, sum(k.diag(A(i,1:i-1))) is the same as summing A(i,j)*k(:,j)
    % for j from 1 to i-1 
    % sum(B,2) returns column vector containing sum of each row of B
    % diag(v) returns matrix with vector v as diagonal.
        k(:, i) = f(t0 + c(i)*h, x0 + h*sum(k(:, 1:i-1)*diag(A(i, 1:i-1)), 2));
    end

    x1 = x0 + h*sum(k*diag(b), 2);
end


%% Class 5 Exercise 3 Functions
% This function plots the stability region for RK4 against explicit Euler
% in the complex plane.

function C5Ex3_plotStability()

% Set up complex plane with mesh of fineness 0.01 squares.
% z = \lambda h (\lambda = A in 1d case)
[x,y] = meshgrid(-3:0.01:1,-3:0.01:3);
z = x+1i*y; 

% RK4 takes the form of a Taylor expansion of order 4 (B from Ex2)
R = abs(1+z+0.5* z.^2 + (1/6) *z.^3 + (1/24)* z.^4); 
% Explicit Euler takes the form of a first order expansion, so converges at
% a narrower range of points
S = abs(1+z); 

figure(4);clf;hold on;
contour(x,y,R,0:0.1:1)
contour(x,y,S,[0.5,1],'r') 

end






