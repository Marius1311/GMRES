%% GMRES - Algorithm
% This is the main skript for the GMRES Algorithm to solve linear Systems of the form Ax = b.

% We need a rectangular, nonsingular Matrix A and a vector b of
% appropriate size. We will try out a couple of different scenarious here:

%% Random Matrix

A = rand(100, 100);
b = rand(100, 1);
x0 = zeros(100, 1);

[iter, resvec, x] = GMRES (A, b, x0);

semilogy(resvec), xlabel('k'), ylabel('log(r_k)');

%% Matrix with few distinct eigenvalues:


X=randn(9,9); 
X=orth(X);
A=X*diag([1,1,4,3,3,4,4,4,3])*X';
b = rand(9, 1);
x0 = zeros(9, 1);

[iter, resvec, x] = GMRES (A, b, x0, 9, 1.e-6);

semilogy(resvec), xlabel('k'), ylabel('log(r_k)');

%% Matrix with given condtion number

% Wel  conditioned if invkappa just less than one
invkappa = 0.9;
A=sprandsym(200,0.1,invkappa,1);
b=ones(200,1);
x0 = zeros(200, 1);

[iter, resvec, x] = GMRES (A, b, x0, 200, 1.e-10);

semilogy(resvec), xlabel('k'), ylabel('log(r_k)');







    


