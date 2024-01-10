function [x,S,r] = LMFsolve(FUN,x)
% MATLABSHARED.TRACKING.INTERNAL.SCENARIO.LMFsolve Solve nonlinear equations
% Xf = matlabshared.tracking.internal.scenario.LMFsolve(FUN,Xo) returns a
% solution, Xf, via a variant of the Levenberg-Marquardt-Fletcher algorithm
% for minimization of a sum of squares of residuals.  FUN is a function
% handle that evaluates the residuals and the corresponding Jacobian
% matrix.  Xo is a column vector of the initial solution.
%
%
%   This function is for internal use only and may be removed in a later
%   release.

%   Copyright 2019 The MathWorks, Inc.
%
%   Adapted from:
%      MathWorks File Exchange LMFsolve version 1.2 
%      Copyright (c) 2007, Miroslav Balda
%      Distributed under BSD License
%
%   Reference:
%   Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least
%   Squares. Rpt. AERE-R 6799, Harwell

%#codegen

n = length(x);
epsx = 1e-7;     % tolerance on solution
epsr = 1e-7;     % tolerance on worst residual error
maxiter = 100*n; % maximum permitted number of iterations

[r,Jtri]= FUN(x);
% J = sparse([1:n 1:n 2:n 1],[2:n 1 1:n 1:n],Jtri,n,n);
% v = J.'*r;

v = mulJt(Jtri,r);
S = r'*r;

% start with a newton step
lambda=0;
lambdac=.75;

iter = 0;

% provide hint to MATLAB Coder
epsx = repmat(epsx,n,1);
epsr = repmat(epsr,n,1);

dx = epsx;

while iter < maxiter && any(abs(dx) >= epsx) && any(abs(r) >= epsr)
    % bump iteration count
    iter = iter+1;
    
    % get new candidate solution
    dx = solveDampenedHessian(Jtri,lambda,v);
    xnew = x-dx;
    
    % evaluate new solution's residual, jacobian, and sum-squared error
    [rnew,Jtrinew] = FUN(xnew);
    Snew = rnew.'*rnew;

    % adjust relaxation coefficient
    [lambda,lambdac] = fletcher(S,Snew,dx,v,Jtri,lambda,lambdac);

    % accept candidate solution when error improves
    if Snew<S
        S = Snew; 
        x = xnew; 
        r = rnew;
        Jtri = Jtrinew;
        v = mulJt(Jtri,r);
    end
end


function [lambda,lambdac] = fletcher(S,Snew,dx,v,Jtri,lambda,lambdac)
% take ratio of actual to predicted reduction in squared error
Adx = mulJtJ(Jtri, dx);
R  = (S-Snew)/(dx.'*(2*v-Adx));

% recompute lambda when ratio is out of tolerance
Rlo = 0.25; 
Rhi = 0.75;

if R>Rhi
    % halve lambda and reset to newton method when lambda is below cutoff
    lambda = lambda/2;
    if lambda<lambdac
        lambda=0;
    end
elseif R<Rlo 
    % compute scaling factor for lambda
    nu = (Snew-S)/(dx.'*v)+2;
    
    % clamp between 2 and 10
    nu = min(max(2,nu),10);
    
    if lambda==0
        % in newton step
        % recompute new values for lambda and cutoff
        % compute via smallest eigenvalue of J'J.
        lambdac = leastEigvJtJ(Jtri,min(40,size(Jtri,1)));
        lambda = lambdac;
        nu = nu/2;
    end
    
    % rescale lambda
    % provide hint to MATLAB Coder
    tmp = nu*lambda;
    lambda = tmp(1);
end


function dx = solveDampenedHessian(Jtri,lambda,v)
n = length(v);

% tridiagonal matrix
a = Jtri(:,1);  % lower  diagonal
b = Jtri(:,2);  % middle diagonal
c = Jtri(:,3);  % upper  diagonal

% symmetric pentadiagonal matrix J'*J
d = (a.^2 + b.^2 + c([n 1:n-1]).^2)*(1+lambda); % diag(J'*J,0)
e = a.*b([2:n 1]) + b.*c;                       % diag(J'*J,1)
f = a.*c([2:n 1]);                              % diag(J'*J,2)

if n<5
    % matrix not pentadiagonal - compute explicitly
    J = accumarray([  1:n    1:n  2:n 1;    ...
                     2:n 1   1:n   1:n  ]', ...
                   [   a;     b;   c    ],[n n]);    
    A = J.'*J;
    dx = (A + diag(lambda*diag(A))) \ v;
elseif a(end) == 0 && c(end) == 0
    dx = sympentdisolve(d, e, f, v);
else
    dx = cycsympentdisolve(d, e, f, v);
end

function [dx,dy] = solveHessian(Jtri,v)
%   solve J'J dx = v.
%   where J is a (cyclic) tridiagonal matrix formed
%   whose diagonals are contained in an n x 3 matrix, Jtri.
%   J = sparse([1:n 1:n 2:n 1],[2:n 1 1:n 1:n],Jtri,n,n)
lower = Jtri(:,1);
center = Jtri(:,2);
upper = Jtri(:,3);

if lower(end) == 0 && upper(end) == 0
    % solve J'*dy = v
    dy = tridisolve(upper, center, lower, v);
    
    % solve J*dx = dy
    dx = tridisolve(lower, center, upper, dy);
else
    % solve J'*dy = v
    dy = cyctridisolve(upper, center, lower, v);
    
    % solve J*dx = dy
    dx = cyctridisolve(lower, center, upper, dy);
end

function x = cyctridisolve(a, b, c, d)
%  Solve the  n x n  tridiagonal system for x:
%
%  [ b(1)  c(1)                             a(n) ] [  x(1)  ]   [  d(1)  ]
%  [ a(1)  b(2)  c(2)                            ] [  x(2)  ]   [  d(2)  ]
%  [       a(2)  b(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        a(n-2) b(n-1) c(n-1) ] [ x(n-1) ]   [ d(n-1) ]
%  [ c(n)                          a(n-1)  b(n)  ] [  x(n)  ]   [  d(n)  ]
%
%  x must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that a(1) and c(n) are not used)

%   Reference:
%      Press, W.H., Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P.
%      (1986) Numerical Recipes, 1st ed., Cambridge University Press.

n = length(d);

gamma = -b(1); % gamma can be arbitrary; -b1 is recommended
b(1) = b(1) - gamma;
b(n) = b(n) - a(n) * c(n)/gamma;

x = tridisolve(a, b, c, d);

u = zeros(n,1);
u(1) = gamma;
u(n) = c(n);             
u(2:n-1) = 0;            

z = tridisolve(a, b, c, u);

fact = (x(1) + a(n)*x(n)/gamma)/(1.0 + z(1) + a(n) * z(n)/gamma);

x = x - fact .* z;


function x = tridisolve(a,b,c,d)
%TRIDISOLVE  Solve tridiagonal system of equations.
%   x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%   [ b(1)  c(1)                                  ] [  x(1)  ]   [  d(1)  ]
%   [ a(1)  b(2)  c(2)                            ] [  x(2)  ]   [  d(2)  ]
%   [       a(2)  b(3)  c(3)                      ] [        ]   [        ]
%   [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%   [                    ...    ...    ...        ] [        ]   [        ]
%   [                        a(n-2) b(n-1) c(n-1) ] [ x(n-1) ]   [ d(n-1) ]
%   [                               a(n-1)  b(n)  ] [  x(n)  ]   [  d(n)  ]
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

%   Copyright 2014 Cleve Moler
%   Copyright 2014 The MathWorks, Inc.

x = d;
n = length(x);

for j = 1:n-1
   mu = a(j)/b(j);
   b(j+1) = b(j+1) - mu*c(j);
   x(j+1) = x(j+1) - mu*x(j);
end

x(n) = x(n)/b(n);
for j = n-1:-1:1
   x(j) = (x(j)-c(j)*x(j+1))/b(j);
end


function x = cycsympentdisolve(d,e,f,b)
%CYCSYMPENTDISOLVE Solve cyclic symmetric pentadiagonal system of equations.
% X = CYCSYMPENTDISOLVE(D,E,F,B) solves A\B for each column of B.
%   D, E, F are the diagonal, off diagonal, and outer diagonal elements
%   of matrix A.  E(end) contains the corner elements, and F(end-1:end)
%   contains the inner corner elements.

%   Reference:
%      Press, W.H., Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P.
%      (1986) Numerical Recipes, 1st ed., Cambridge University Press.

n = length(b);

% compute U
u = zeros(n,4);
u(1,1) = 1;
u(2,2) = 1;
u(n-1,3) = 1;
u(n,4) = 1;

% sympentdisolve ignores the cyclic entries
zy = sympentdisolve(d,e,f,[u b]);

% extract z and y from solver
z = zy(:,1:4);
y = zy(:,5);

% compute V
v = zeros(n,4);
v(n-1,1) = f(n-1); %A(1,end-1);
v(n,1) = e(n);     %A(1,end);
v(n,2) = f(n);     %A(2,end);
v(1,3) = f(n-1);   %A(end-1,1);
v(1,4) = e(n);     %A(end,1);
v(2,4) = f(n);     %A(end,2);

% Eqns 2.7.19, 2.7.21
x = y - z * ((eye(4) + v'*z) \ (v'*y));


function x=sympentdisolve(d,e,f,b)
%SYMPENTDISOLVE Solve symmetric pentadiagonal system of equations.
% X = SYMPENTDISOLVE(D,E,F,B) solves A\B for each column of B.
%   D, E, F are the diagonal, off diagonal, and outer diagonal elements
%   of matrix A.

%   Copyright 2019 The MathWorks, Inc.
%
%   Adapted from:
%      MathWorks File Exchange Fast Pentadiagonal System Solver
%      Copyright (c) 2009, Greg von Winckel
%      Distributed under BSD License

%   Reference:
%      G. Engeln-Muellges, F. Uhlig, "Numerical Algorithms with C"
%      Chapter 4. Springer-Verlag Berlin (1996)

[N,P] = size(b);
x = zeros(N,P);
    
alpha = zeros(N,1);
gamma = zeros(N-1,1);
delta = zeros(N-2,1);

% Factor A=LDL'
alpha(1) = d(1);
gamma(1) = e(1)/alpha(1);
delta(1) = f(1)/alpha(1);

alpha(2) = d(2)-e(1)*gamma(1);
gamma(2) = (e(2)-f(1)*gamma(1))/alpha(2);
delta(2) = f(2)/alpha(2);

for k=3:N-2
    alpha(k) = d(k)-f(k-2)*delta(k-2)-alpha(k-1)*gamma(k-1)^2;
    gamma(k) = (e(k)-f(k-1)*gamma(k-1))/alpha(k);
    delta(k) = f(k)/alpha(k);
end

alpha(N-1) = d(N-1)-f(N-3)*delta(N-3)-alpha(N-2)*gamma(N-2)^2;
gamma(N-1) = (e(N-1)-f(N-2)*gamma(N-2))/alpha(N-1);
alpha(N)   = d(N)-f(N-2)*delta(N-2)-alpha(N-1)*gamma(N-1)^2;

% Update Lx=b, Dc=z
z=zeros(N,P);
z(1,:) = b(1,:);
z(2,:) = b(2,:)-gamma(1)*z(1,:);

for k=3:N
    z(k,:) = b(k,:) - gamma(k-1)*z(k-1,:) - delta(k-2)*z(k-2,:);
end

c = bsxfun(@rdivide, z, alpha);

% Backsubstitution L'x=c
x(N,:) = c(N,:);
x(N-1,:) = c(N-1,:) - gamma(N-1)*x(N,:);

for k=N-2:-1:1
    x(k,:) = c(k,:) - gamma(k)*x(k+1,:) - delta(k)*x(k+2,:);
end

function lev = leastEigvJtJ(Jtri,n)
%LEASTEIGJTJ - solve for least eigenvalue of J'J
%  LEV = LEASTEIGJTJ(Jtri,N) approximates the least eigenvalue of J'J via
%  N Lanczos iterations.  J is a (cyclic) tridiagonal matrix formed
%  whose diagonals are contained in an n x 3 matrix, Jtri.
%  J = sparse([1:n 1:n 2:n 1],[2:n 1 1:n 1:n],Jtri,n,n)

bet = 0;
m = size(Jtri,1);

oldq = zeros(m,1);
b = ones(m,1);

q = b / norm(b);
alpha = zeros(n,1);
beta = zeros(n,1);

for i=1:n
    v = solveHessian(Jtri,q); 
    alpha(i) = q'*v;
    v = v - bet*oldq - alpha(i)*q;
    bet = norm(v);
    oldq = q;
    q = v / bet;
    beta(i) = bet;
end

% Create symmetric tridiagonal matrix
% use full matrix instead of sparse for MATLAB Coder.
T = accumarray(    [1:n-1      1:n         2:n;
                    2:n        1:n       1:n-1]', ...
                [beta(1:n-1); alpha; beta(1:n-1)],[n n]);
lev = 1./min(abs(eig(T)));


function y = mulJt(Jtri, x)
% compute y = J*x
% where J = sparse([1:n 1:n 2:n 1],[2:n 1 1:n 1:n],Jtri,n,n)
n = length(x);
lower = Jtri(:,1);
center = Jtri(:,2);
upper = Jtri(:,3);

y = center.*x + lower.*x([2:n 1]) + upper([n 1:n-1]).*x([n 1:n-1]);

function y = mulJtJ(Jtri, x)
% compute y = J'*J*x
% where J = sparse([1:n 1:n 2:n 1],[2:n 1 1:n 1:n],Jtri,n,n)
n = length(x);
lower = Jtri(:,1);
center = Jtri(:,2);
upper = Jtri(:,3);

z = center.*x + upper.*x([2:n 1]) + lower([n 1:n-1]).*x([n 1:n-1]);
y = center.*z + lower.*z([2:n 1]) + upper([n 1:n-1]).*z([n 1:n-1]);
