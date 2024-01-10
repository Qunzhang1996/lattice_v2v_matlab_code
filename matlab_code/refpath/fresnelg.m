function z = fresnelg(x,dk,k,theta)
%MATLABSHARED.TRACKING.INTERNAL.SCENARIO.FRESNELG Generalized Fresnel integral.
%
%   This function is for internal use only and may be removed in a later
%   release.
%
%   Z = MATLABSHARED.TRACKING.INTERNAL.SCENARIO.FRESNELG(L, DK, K, THETA)
%   returns the complex Fresnel integral:
%
%              /L
%              |
%        Z  =  |  exp(1i * ((DK/2)*s.^2 + K*s + THETA)) ds
%              |
%              /0
%    
%   where DK, K, and THETA, are real scalar constants, and L is real
%   N-dimensional array of points over which to evaluate the integral.
%
%   If omitted, DK = pi, K = 0, THETA = 0.

%#codegen

%   Copyright 2017-2021 The MathWorks, Inc.

% coder.noImplicitExpansionInFunction;

if nargin<2
    dk = pi;
end

if nargin<3
    k = 0;
end

if nargin<4
    theta = 0;
end

validateattributes(x,{'double'},{'real'},'fresnelg','x',1);
validateattributes(dk,{'double'},{'real','scalar','finite'},'fresnelg','dk',2);
validateattributes(k,{'double'},{'real','scalar','finite'},'fresnelg','k',3);
validateattributes(theta,{'double'},{'real','scalar','finite'},'fresnelg','theta',4);

if dk ./ (k.*k) > 1e-6 %#ok<BDSCA>
    z1 = fresnel(sqrt(dk./pi)*x+k./sqrt(pi*dk));
    z0 = fresnel(k./sqrt(pi*dk));
    z = sqrt(pi./dk)*exp(1i*(theta-k.*k./(2*dk)))*(z1-z0);
elseif dk ./ (k.*k) < -1e-6 %#ok<BDSCA>
    z1 = fresnel(sqrt(-dk./pi)*x-k./sqrt(-pi*dk));
    z0 = fresnel(-k./sqrt(-pi*dk));
    z = conj(sqrt(-pi./dk)*exp(-1i*(theta-k.*k./(2.*dk)))*(z1-z0));
else
    z = fresnelgsma(x,dk,k,theta);
end

% overwrite nearly straight segments
% Implicit expansion disabled here
izero = find(abs(dk).*x.^2 < 1e-3 & abs(k.*x) < 1e-3);
if ~isempty(izero)
    z(izero) = fresnelgzero(x(izero),dk,k,theta);
end

function z = fresnelgsma(x,dk,k,theta)
% compute for large or unstable values of k^2/dk
if k==0 % then dk==0
    z = x .* exp(1i .* theta);
else
    z = fresnelgsmaexp(x,dk,k,theta);
end

function z = fresnelgsmaexp(x,dk,k,theta)
coder.noImplicitExpansionInFunction;

C = 0.5i*dk ./ k.^2;

nikx = -1i*k .* x;
e = exp(-nikx);
t = 1 - e;
m = -e;
c = ones(size(C),'like',C);

% provide hint to MATLAB coder for size of s.
s = c .* t;

for n=1:10
    % Implicit expansion disabled here
    m = m .* nikx;
    r = 2*n .* m;
    % Implicit expansion disabled here
    m = m .* nikx;
    % Implicit expansion disabled here
    t = (2*n-1)*(2*n) .* t + r + m;
    c = c .* C./-n;
    % Implicit expansion disabled here
    s = s + c .* t;
end
    
z = s .* 1i./k .* exp(1i*theta);

function z = fresnelgzero(x,dk,k,theta)
% coder.noImplicitExpansionInFunction;

N = 5;

% expand exp(i*(dk/2).*x^2)
% assume A(:,0) = 1.
% A(:,k) = (1i*(dk/2)*x.^2))^k/k!
a = 0.5i*dk*x(:).^2;
A = repmat(a,1,N);
A = bsxfun(@rdivide,A,1:N);
A = cumprod(A,2);

% expand exp(i*k.*x)
% assume B(:,0) = 1.
% B(:,k) = (1i*k.*x)^k/k!
b = 1i*k.*x(:);
B = repmat(b,1,N);
B = bsxfun(@rdivide,B,1:N);
B = cumprod(B,2);

z = ones(length(x),1,'like',complex(x(1)));
for i=1:N
    % Implicit expansion disabled here
    z = z + A(:,i)./(2*i+1);
end

for j=1:N
    % Implicit expansion disabled here
    z = z + B(:,j)./(j+1);
end

% do cross terms
for i=1:N-1
    for j=1:N-2*i
        % Implicit expansion disabled here
        z = z + A(:,i).*B(:,j)./(2*i+j+1);
    end
end

% Implicit expansion disabled here
z = z .* x(:) .* exp(1i*theta);
