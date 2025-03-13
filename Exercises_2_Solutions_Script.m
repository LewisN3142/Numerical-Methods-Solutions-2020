% Displays numeric values of type 'double' to 15 decimal places and 
% those of type 'single' to 7 decimal places.
format long;

%% Class 2 Exercise 1 Script
% This runs the C2Ex1_LUSolver function on a randon 4x4 matrix and 4x1
% column vector, returning an approximate solution of Ax = b, as well as the
% residual.

A_Example = rand(4)
b_Example = rand(4,1) 
[x_Example,r_Example] = C2Ex1_LUSolver(A_Example,b_Example)


%% Class 2 Exercise 2 Script
% This runs the C2Ex2_hilbertProblem function for eps=0 and eps=10^-3 for
% the 5x5 Hilbert matrix and plots the runtime of the perturbed and
% unperturbed Hilbert problems of sizes 10:10:1000.

[xa,ta] = C2Ex2_hilbertProblem(5,0)
[xb,tb] = C2Ex2_hilbertProblem(5,1e-3)
c = cond(C2Ex2_hilbertMatrix(5))

C2Ex2_hilbertPlot(10,10,100,1e-3)

%% Class 2 Exercise 3
% Maths exercise, see pdf for solution.


%% Class 2 Exercise 1 Functions
% A series of functions to numerically compute algebraic systems of
% equations of the form Ax = b.

% This function solves the equation Ax=b 
% LU decomposition gives LUx = b
% forward substitution gives Ux = y
% backward substitution gives x 
% We also return r, the residual of Ax = b, to see how accurate the
% approximate solution is.

function [x,r] = C2Ex1_LUSolver(A,b)

% Begin by checking length of vector and dimension of matrix match
if size(A,1) ~= size(b)
    error('Dimensions of matrix and vector not compatible')
end

[L,U] = C2Ex1_LUDecomp(A);
y = C2Ex1_fsub(L,b);
x = C2Ex1_bsub(U,y);
r=max(abs(A*x-b));
end


% This function performs the LU decomposition process on a matrix A,
% returning two triangular matrices L and U.

function [L,U] = C2Ex1_LUDecomp(A) 

% Begin by checking if input argument is square matrix.
s=size(A); 
if s(1) ~= s(2)
   error('Not a square matrix') 
end

n=s(1);

% One could alternatively use nested for loops here, but this uses matlab
% efficiency.
for k = 1:n-1 
    % Set rho to be a vector of values k+1,k+1,... n
    rho = k+1:n; 
    % Update last n-k values in kth column
    A(rho, k) = A(rho, k)/A(k, k);
    % This re-evaluates the bottom right n-k by n-k submatrix
    A(rho, rho) = A(rho, rho) - A(rho, k)*A(k, rho); 
end
    
% Extract lower triangular part of new A and set diagonal to 1s
L = tril(A, -1) + eye(n); 
% Extract upper triangular part (inc. diagonal)
U = triu(A);

end


% This function performs forward substitution on a lower triangular matrix
% A and column vector b, returning a vector x which solves Ax = b.
% As our arguments come from LU decomposition, we know they are of the correct
% type -- how would you check this if it wasn't guaranteed?

function x = C2Ex1_fsub(A,b) 
n=length(b); 
x = zeros(n,1);
s = zeros(n,1);
for i = 1:n 
    x(i)=0; 
    s(i)=0;
    for j=1:i-1 %
        s(i)=s(i)+A(i,j)*x(j); 
    end
    x(i)=(b(i)-s(i))/A(i,i); 
end
end


% This function performs backward substitution on a upper triangular matrix
% A and column vector b, returning a vector x which solves Ax = b.
% As our arguments come from LU decomposition, we know they are of the correct
% type -- how would you check this if it wasn't guaranteed?

function x = C2Ex1_bsub(A,b) 
n=length(b);
x=zeros(n,1); 
for i = n:-1:1 
    s(i)=0;
    for j=i+1:n
        s(i)=s(i)+A(i,j)*x(j); 
    end
    x(i)=(b(i)-s(i))/A(i,i);
end
end


%% Class 2 Exercise 2 Functions
% A series of functions to solve the perturbed and unperturbed Hilbert
% problem, then plot their computation times.

% This function takes a range of matrix sizes s:step:f and a perturbation,
% and plots the runtime for the perturbed Hilbert problem against the
% unperturbed version.

function [N,T1,T2] = C2Ex2_hilbertPlot(s,step,f,eps)
% Set the range of matrix sizes to test
N = s:step:f;
lenN = length(N);

% Initialise vectors to store the runtimes
T1 = zeros(length(N),1);
T2 = zeros(length(N),1);


for i = 1:lenN
    % Capture the times it takes in T1 and T2
    [~,T1(i)] = C2Ex2_hilbertProblem(N(i),0);
    [~,T2(i)] = C2Ex2_hilbertProblem(N(i),eps);
end

figure(1); %Straightforward plot of both sets of timestamps
plot(N,T1,N,T2);
title('Runtime')
legend('Standard','Perturbed')
xlabel('Matrix Size')
ylabel('Runtime (s)')

figure(2); %Plot graph with logarithms of both scales
loglog(N,T1,N,T2);
title('Loglog Runtime')
legend('Standard','Perturbed')
xlabel('Matrix Size')
ylabel('Runtime (s)')

end


% This function takes a number n (the size of a square matrix) and
% perturbation eps and outputs the approximate solution x to the Hilbert
% problem perturbed by epsilon, as well as the time t it takes to solve.

function [x,t] = C2Ex2_hilbertProblem(n,eps)
% Generate the Hilbert matrix of size n
H = C2Ex2_hilbertMatrix(n); 
% Generate the vector b (note this is derived from H)
b = H*ones(n,1);
if eps ~= 0
b(1) = b(1)+eps;
end
% Start timing and solve Hx=b
tic;
x = C2Ex1_LUSolver(H,b);
t = toc;
end


% This function makes the n by n Hilbert matrix with entries H(i,j)

function H = C2Ex2_hilbertMatrix(n)
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            H(i, j) = 1/(i + j - 1);
        end
    end
end