function [iter, resvec, x] = GMRES (A, b, x0, varargin)

% We need to check whether A is singular and whether the dimensions are
% appropriate: (One shouldn't do this in reality, it is much much too expensive)

%if det(A) == 0
%    disp('A singular')
%end

[n, n] = size(A);
if n ~= length(b) | n ~= length(x0)
    disp('Wrong dimensions')
end

%% Set parameter values

% We need to set parameters such at the maximum number of iterations and
% the tolerance level

if nargin == 3
    TOL = 1.e-3;
    iterMax = n;
elseif nargin == 4
        iterMax = varargin{1};
        TOL = 1.e-4;
elseif nargin == 5
        iterMax = varargin{1};
        TOL = varargin{2};
end

%% Initialise values

% Create r0
r0 = b - A*x0;
% Create first basis vector for Krylov subspace
V = zeros(n, n);
V(:, 1) = r0/norm(r0, 2);
% Set loop counter
k = 1;
% Create the Matrix H
H = zeros(2, 1);
% Matrices for QR factorisation
R = H;
G = cell(n, 1);
resvec = zeros(n, 1);
resvec(1, 1) = norm(r0, 2);
RES = resvec(1, 1);

%% Main routine

while RES > TOL && k <= iterMax
 
    %---------------------------
    % Do step k of the Arnolid Algorithm
    w = A*V(:, k);
    h = zeros(k+1, 1);
    for j = 1:k
        h(j, 1) = V(:, j)'*w;
        w = w - h(j, 1)*V(:, j);
    end
    h(k+1, 1) = norm(w, 2);
    
    % if we are at k = n, then we have already found a basis for the whole
    % R^n, so we do not need any more voctors in V
    if k ~= n
        V(:, k+1) = w/h(k+1, 1);
    end
    %---------------------------
    
    %---------------------------
    % Solve linear least squares problem to find residual
    if k == 1
        H = h;
    else
        H = [[H; zeros(1, k-1)], h];
    end
    
    % We need QR fact. of H. Use givens rotations.
    %Store H in R
    R = H;
    
    for j = 1:k
        % We only actually have to compute a givens rotation in the last
        % case...
        if j == k
           G{k, 1} = eye(k+1); 
           [c,s] = givensrotation( R(k,k),R(k+1,k) );
           G{k, 1}([k, k+1],[k, k+1]) = [c -s; s c]; 
           Givens = G{k, 1};
        % Otherwise, we use one of the olds ones and adkust the dimensions
        else
            Rest = eye(k+1);
            Givens = [G{j, 1}, Rest(1:(j+1), (j+2):(k+1)); Rest(j+2:k+1, 1:k+1)];
        end
        % We update R
        R = Givens'*R;
        % And we update Q
        if j == 1
            Q = Givens;
        else
            Q = Q*Givens;
        end
    end
    
    %We will now reduce the least squares problem to a problem with unique
    %solution:
    Reff = R(1:k, :);
    c = norm(r0, 2) * Q' * eye(k+1, 1);
    ceff = c(1:k);
    % This gives me the linear least squares error, which is my residual
    resvec(k+1, 1) = abs(c(k+1));
    RES = resvec(k+1, 1);
    % We have not yet solved the linera least squares problem, but we do
    % not have to. We will only do this once the residual is small enough.
    
    % update k
    k = k+1;
end


%% Compute solution

% If we are happy with teh size of the residual, we can actually compute
% the solution of the linear least squares problem and therewith the
% solution of the problem:

y = zeros(k-1, 1);

% We will simply use back-substitution:
for i = (k-1):-1:1
         back = 0;
    for j = (k-1):-1:(i+1)
         back = back + Reff(i, j) * y(j, 1);
    end
    y(i, 1) = (ceff(i) - back)/Reff(i, i);
end

% We have to adjust the size of V:
Veff = V(:, 1:(k-1));

% Now, the solution to the problem is given by:
x = x0 + Veff*y;

% We shorten the residual vector:
resvec = resvec(1:k);

iter = k-1;

end

