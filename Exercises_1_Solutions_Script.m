% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long;

%% Class 1 Exercise 1 Script
% This runs the C1Ex1_CubicRoot function on several initial
% conditions and returns the three roots, as well as the number of
% iterations required to reach the desired accuracy.

[Q1i_root1,Q1i_root2,Q1i_root3,Q1i_iterations] = C1Ex1_CubicRoot(3.375, 10^-6)
[Q1ii_root1,Q1ii_root2,Q1ii_root3,Q1ii_iterations] = C1Ex1_CubicRoot(8, 10^-6)
[Q1iii_root1,Q1iii_root2,Q1iii_root3,Q1iii_iterations] = C1Ex1_CubicRoot(1331, 10^-6)


%% Class 1 Exercise 2 Script
% This runs the C1Ex2_2Norm function on a random vector 'x' of random
% length 'length' and outputs the 2-norm as a real number 'y.' For comparison,
% we also compute the 2-norm using the inbuilt functions.

length = randi([1 30],1,1)
x = rand(length, 1)
our2Norm = C1Ex2_2Norm(x)
matlab2Norm = norm(x,2)

%% Class 1 Exercise 3 Script
% Arguments for the mystery function
A = [2, -1, 0; -1, 2, -1; 0, -1, 2]
x = [1; 1; 1]

C1Ex3_MysteryFunction(x,A,10)

% Compare with the eigenvalues of A -- the algorithm approximates the largest.
[V,D] = eig(A)


%% Class 1 Exercise 4 Script

% Compute the determinant of a random matrix using the Laplace expansion
% and the inbuilt Matlab function
A_Example = rand(4)
det(A_Example)
C1Ex4_detLaplace(A_Example)

% Plot the runtime of the Laplace expansion against the inbuilt function
% for square matrices of size 1 to 10.
[t1,t2]=C1Ex4_detTime(10)


%% Class 1 Exercise 1 Function
% This function computes the cubic root using the second algorithm 
% It takes a number 'a' and a desired accuracy 'tolerance' and returns all
% three cubic roots, as well as the number of iterations to reach the
% desured accuracy.

function [x1,x2,x3,iteration] = C1Ex1_CubicRoot(a,tolerance)

% Initialise local variables.
t = tolerance;
x1 = a;
iteration = 0;

% Second algorithm for computing solutions of a cubic.
while(abs(x1^3-a)>= t)
    x1=(2*x1+a/x1^2)/3;
    iteration = iteration+1;

    % Uncomment the line below to print each iteration to the command window
    %fprintf('Iteration %d: x = %f\n', iteration, x1);
end
% Compute the other solutions using roots of unity.
x2=x1*exp(1i*2*pi/3);
x3=x1*exp(1i*4*pi/3);
end

%% Class 1 Exercise 2 Function
% This function takes a vector 'x' and outputs its Euclidean norm (2-norm)
% as 'norm.'

function norm = C1Ex2_2Norm(x)

% Here we have a first example of argument checking and future-proofing. We
% make sure to check whether or not the argument is a vector and if it is,
% whether or not we have a row or column. It does not matter here, but in
% some scripts, the orientation of the vector arguments will be important.
sizeX = size(x);
if sizeX(1)==1
    fprintf('You have entered a row vector\n')
    fprintf('The length of your vector is %d \n',sizeX(2))
    n = sizeX(2);
elseif sizeX(2)==1
    fprintf('You have entered a column vector\n')
    fprintf('The length of your vector is %d \n',sizeX(1))
    n= sizeX(1);
else
    error('You have entered an argument which is not a vector \n')
end
norm = 0;

for i = 1:n
   norm=norm+x(i)^2; 
end
norm=sqrt(norm);
end

%% Class 1 Exercise 3 Function
% This function takes a matrix 'A' and a vector 'x' and performs our
% mystery algorithm for 'iter' iterations, before outputting some number
% 'y'

function y = C1Ex3_MysteryFunction(x,A,iter)

a = 0;
for i = 1:iter
   x = A*x;
   x = x*1/norm(x);
   a = (dot(A*x,x))/(dot(x,x));
end
y=a;
end


%% Class 1 Exercise 4 Function 
% A pair of functions to compute the determinant of a matrix using the
% Laplace expansion and compare the runtime of this operation to the
% inbuilt matlab function.

% This function takes a number s and computes the determinant of a random
% square matrix of size 1 up to s, using both the Laplace expansion and the
% inbuilt matlab function, then plots the runtime for these functions.

function [t1,t2] = C1Ex4_detTime(s) 
t1 = zeros(s,1); 
t2 = zeros(s,1);

% Loop through the selected range of matrix sizes
for i = 1:s
     A = rand(i);
     
     % Compute the determinant using the inbuilt function
     tic;
     d1 = det(A);
     t1(i) = toc;
     
     % Compute the determinant using our custom function
     tic;
     d2 = C1Ex4_detLaplace(A);
     t2(i) = toc;
  
     A;
     fprintf('Matrix size: %x\n In built determinant: d1 = %f, t1 = %f\n Laplace: d2 = %f, t2 = %f\n',i,d1,t1(i),d2,t2(i))    
 end

 figure(1);clf;hold on;
 plot(t1);
 plot(t2);
 title('Plot of runtimes for Laplace and Matlab computing determinant')
legend('MATLAB','Laplace')
xlabel('Size of matrix')
ylabel('Runtime (s)')

end


% This function computes the determinant of a matrix A using the Laplace
% expansion from lectures. 

function s = C1Ex4_detLaplace(A)

% We first check that the matrix is square and therefore has a determinant.
n = size(A);
if n(1) ~= n(2)
    error('Not a square matrix')
end

% The determinant of a scalar is the scalar
if n(1)==1
    s=A;
    return;
end

% The determinant of the matrix A is found by taking the i-th entry in the first row of A, 
% multiplying this by the determinant of the minorformed by removing row 1 
% and column i, then summing over i (up to a sign).
s = 0;
    for i = 1:n
        ind = [1:i-1, i+1:n];
        s = s + (-1)^(i+1) * A(1, i) * C1Ex4_detLaplace(A(2:n, ind));
    end
end
