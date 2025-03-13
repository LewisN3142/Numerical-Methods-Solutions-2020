% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long;

%% Class 4 Exercise 1 Script
% This exercise is purely mathematical. See pdf for solution.


%% Class 4 Exercise 2 Script
% This script executes the functions to approximate f(x) = exp(x)
% on the interval [-1,1] using Gauss' best approximation against the 
% basis {1, x, x^2} and plots the result against f(x) itself.

C4Ex2_bestApprox(@(x) exp(x),-1,1)

%% Class 4 Exercise 3 Script
% Given a function f and its exact integral fInt, we compute Simpsons rule
% for a range of subdivisions of [5,10] (namely 10:10:500) and show the
% rule converges as the number of subdivisions increases, while error
% decreases with subdivision width.

f = @(x) 1/x;
fInt = @(x) log(x);
C4Ex3_simpsonPlot(f,fInt,5,10,10:10:500)



%% Class 4 Exercise 2 Functions
% This function takes an interval [a,b] and returns the coefficients \alpha
% in an expansion against basis \phi, approximating the function f

function alpha = C4Ex2_bestApprox(f,a,b)

% Set basis to approximate against (could also use Fourier basis
% @(x)(ones), @(x) cos(x), @(x) sin(x), @(x) cos(2*x), @(x) sin(2*x)
phi = {@(x) ones(size(x)), @(x) x, @(x) x.^2 };  

% Size of basis we expand against
n = size(phi,2);
% Matrix to store the L^2 inner products of the basis vectors 
C = zeros(n);
% Matrix to store the L^2 inner product of f with basis vectors
r = zeros(n,1);

% Evaluate inner products using integral (quad depreciated, could also use Simpson's rule from exercise 3)
for i=1:n
    for j=1:n
    h = @(x) phi{i}(x).*phi{j}(x);
    C(i,j) = integral(h,a,b);
    end
    h = @(x) phi{i}(x).*f(x); 
    r(i) = integral(h,a,b);
end

% Solve the system of equations C\alpha = r to obtain coefficients in
% expansion f = \alpha_0 + \alpha_1 x + \alpha_2 x^2
alpha = C\r; 

% Constructs approximating polynomial for plot
g = @(x) 0;
for i=1:n
   g = @(x)alpha(i)*phi{i}(x) + g(x); 
end

%Plot approximating polynomial against function
figure; clf; hold on; 
    x=a:0.1:b;
    plot(x,f(x));
    plot(x,g(x));
    legend('Function','Approximation');
end


%% Class 4 Exercise 3 Functions
% A pair of functions to compute the integral of a function f over an
% interval [a,b] using Simpson's rule and compare the result to the exact
% integral.

% This function takes some f, its exact integral fInt, an interval [a,b],
% and a list of numbers of steps N and uses Simpson's rule to compute the
% integral of f with n \in N steps, then compares this to fInt, producing
% plots of the change in error with increase in n.

function I = C4Ex3_simpsonPlot(f,fInt,a,b,N)

% Initialise vectors for Simpson rule output and error
n = length(N);
I = zeros(n,1); 
E = zeros(n,1);
% Carry out Simpson rule I  and error E for each number of steps n in N 
for j=1:n 
   I(j) = C4Ex3_Simpson(f,a,b,N(j));
   E(j) = abs(I(j) - fInt(b) + fInt(a)); 
end

% Plot Simpson's rule against number of intervals, show convergence as n increases
figure(1); 
plot(N,I);
title('Simpson Rule for varying number of intervals.')
xlabel('Number of intervals n')
ylabel('Result of Simpson Rule')


% Plot error against interval width, showing decrease as width h decreases
figure(2); 
loglog((b-a)./N,E);
title('Error of Simpson Rule with varying interval width')
xlabel('Interval width h (log)')
ylabel('absolute error (log)')

end


% This function computes the integral of f on the interval [a,b] using
% Simpson's rule with n intervals.

function I = C4Ex3_Simpson(f,a,b,n)
% Failsafe in case n is not even
if mod(n,2) ~= 0 
    error('Function requires even number of intervals - check last argument even')
end

% Set interval width
h=(b-a)/n; 

% Compute integral
I=0; 
for j=1:2:n-1
    I=I+h/3 * (f(a+(j-1)*h)+4*f(a+j*h)+f(a+(j+1)*h)); 
end
end
