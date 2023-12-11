%% Class 1 Exercise 1 Function
% This function computes the cubic root using the second algorithm 
% It takes a number 'a' and a desired accuracy 'tolerance' and returns all
% three cubic roots, as well as the number of iterations to reach the
% desured accuracy.

function [x1,x2,x3,iteration] = C1Ex1_CubicRoot(a,tolerance)

% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long; 

% Initialise local variables.
t = tolerance;
x1 = a;
iteration = 0;

% Second algorithm for computing solutions of a cubic.
while(abs(x1^3-a)>= t)
    x1=(2*x1+a/x1^2)/3;
    iteration = iteration+1;
end
% Compute the other solutions using roots of unity.
x2=x1*exp(1i*2*pi/3);
x3=x1*exp(1i*4*pi/3);

end

