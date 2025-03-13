% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long;

%% Class 3 Exercise 1 Script
% This runs the C3Ex1_fixedPointPlot and C3Ex1_fixedPoint functions,
% showing that there is a root in [1.75,2], that our scheme is a
% contraction, and that the error decreases logarithmically.

C3Ex1_fixedPointPlot(1.75,2);
C3Ex1_fixedPoint(1.8,20);


%% Class 3 Exercise 2 Script
% Executes C3Ex2_intervalBisection function on the scheme defined by F
% Carries out the first 10 iterations with no required tolerance to stop
% earlier.

F = @(x) x - 1 - 1./x - 1./x.^2;
[a,b] = C3Ex2_intervalBisection(F,1.75,2,0,10) 



%% Class 3 Exercise 3 Script

% This is a plot of the Lagrange interpolating polynomial computed by hand
% in the pdf solutions.
F = @(x) 6.*x.^2./5 - 3.*x./5 - 1;
figure(5);clf;hold on;
fplot(F,[-2,3]);
title('Plot of interpolating polynomial')
xlabel('x')
ylabel('interpolating polynomial p')

% We can also generate the same plot by computing the interpolating
% polynomial numerically. Here we compute the value of the interpolating
% polynomial at values of x in steps of 0.1 from -2 to 3.
figure(6)
plot(-2:0.1:3,lagrange([-2,0,3],[5,-1,8],-2:0.1:3)) 
title('Plot of interpolating polynomial')
xlabel('x')
ylabel('interpolating polynomial p')


%% Class 3 Exercise 4 Script
% Maths exercise, see pdf for solution.


%% Class 3 Exercise 1 Function
% A series of functions to plot the fixed point scheme F and its derivative
% dF, as well as the error between the results of the scheme at step i and
% the value provided by Matlab's inbuilt function.

% This function takes the interval [a,b] in which we wish to find a root of
% x_{k+1} = F(x_k) and highlights whether or not there is a root and if F
% is a contraction on the interval.

function C3Ex1_fixedPointPlot(a,b)

% Symbolic functions for our fixed point scheme function F and its
% derivative dF
F = @(x) 1+1./x +1./x.^2; 
dF = @(x) -1./x.^2 - 2./x.^3;

% We plot the function F and the line y=x.
% Note the curves intersect at one point in the interval. This is the root
% of our fixed point scheme.
figure(1);clf;hold on;
fplot(@(x)x,[a,b]);
fplot(F,[a,b]);
title('Plot of function F')
legend('x','1+1/x+1/x^2')
xlabel('x')
ylabel('y')

% Plot the derivative of F (dF) and show that its modulus is bounded by 1,
% which implies F is a contraction and our scheme converges.
figure(2);clf;hold on; 
fplot(dF,[a,b]);
fplot(1,[a,b]);
fplot(-1,[a,b]);
title('Plot of function dF')
legend('x','-1/x^2 - 2/x^3')
xlabel('x')
ylabel('y')
end


% This function takes an initial guess x for what the root might be
% and a number of iterations of the fixed point scheme to perform.
% The function we perform the scheme for is x^3 - x^2 -x -1 = 0 .

function m = C3Ex1_fixedPoint(xInit,iter) 

% Symbolic functions for our fixed point scheme function F and its
% derivative dF
F = @(x) 1+1./x +1./x.^2; 
dF = @(x) -1./x.^2 - 2./x.^3; 

% Gets roots of our polynomial x^3 - x^2 - x - 1 = 0
s = roots([1,-1,-1,-1]);
% Root of polynomial we are looking for (given by inbuilt function) 
x_star = s(1)

% Vector to store result of each iteration
x = zeros(iter,1);
x(1) = xInit;
fprintf('x_0 = %f\n',x(1));
% Error -- how close our scheme was to Matlab's
e = zeros(iter,1);
e(1) = abs(xInit - x_star);
fprintf('e_0 = %f\n',e(1));

% Compute iter iterations of fixed point scheme
for i = 1:iter 
    x(i+1) = F(x(i));
    fprintf('x_%d = %f\n',i, x(end));
    e(i+1) = abs(x(i) - x_star);
    fprintf('e_%d = %f\n',i, e(end));
end

%plots e against i, with e in a log scale (error deceases exponentially)
figure(3);
semilogy(e); 
title('Plot of error of scheme')
xlabel('number of iterations')
ylabel('error (log)')

%plots the pairwise ratio of errors to check if we have convergence (does converge)
figure(4);
plot(abs(e(2:end)./e(1:end-1))); 
title('Plot of ratios of errors at each step')
xlabel('number of iterations')
ylabel('step i+1 error / step i error')

%Find value of derivative at fixed point
m = (dF(x_star))
end

%% Class 3 Exercise 2 Function
% A function which executes the bisection method for a given function f in
% the interval [a,b], with maximum number of iterations iter and tolerance
% e. Note f should be a symbolic function.

function [a,b]= C3Ex2_intervalBisection(f,a,b,e,iter)  

% Check whether iteration is necessary (error if root at given endpoint)
% and possible (error if two roots in interval provided)
if f(a)==0  
    fprintf('Root is at the endpoint %f\n',a);
    error ('Iteration not needed');
end
if f(b)==0  
    fprintf('Root is at the endpoint %f\n',b);
    error ('Iteration not needed');
end
if (f(a)*f(b)>0)
    error ('Even number of roots in interval provided');
end

% Count number of iterations and stop bisection if equal set limit
counter=0;
% Repeat bisection until tolerance or iteration cap met.
while ((b-a)>=e)
    c=(a+b)/2;
    
    %If the midpoint is a root, print that the root was found and exit while loop.
    if f(c)==0  
        fprintf('After %d iterations, an exact root was found at %f\n', counter,c);
        break;
    end
    if (f(a)*f(c)<0) %If the root is in the left half
        b=c; 
    else %If the root is in the right half
        a=c;
    end
    counter=counter+1;
    fprintf('Iteration number %d, root in interval: [%f,%f]\n',counter,a,b);
    
    if counter==iter
        break;
    end
end
end


%% Class 3 Exercise 3 Function
% This function computes the values of the Lagrange polynomials and
% interpolating polynomial at each of the points in the range x = -2:0.1:3
% using the vectors pointx and pointy from the provided table.

function y=lagrange(pointx,pointy,x)

% Store the length of the vector of provided values
n=size(pointx,2);
% Initialise a matrix to store the computed values of the Lagrange polynomials
% at each point in x
L=ones(n,size(x,2)); 

% Check that the table of provided points is complete
if (n~=size(pointy,2))
   error('x and y must be vectors of same length');
else
   y=0;
   for i=1:n
       % Compute the Lagrange polynomials at each point in x
      for j=1:n
         if (i~=j)
            L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
         end
      end
      % Compute the interpolating polynomial at each point in x
      y=y+pointy(i)*L(i,:);
   end
end
end
