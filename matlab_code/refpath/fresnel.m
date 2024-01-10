function z = fresnel(x)
%MATLABSHARED.TRACKING.INTERNAL.SCENARIO.FRESNEL Fresnel integral.
%
%   This function is for internal use only and may be removed in a later
%   release.
%
%   Z = MATLABSHARED.TRACKING.INTERNAL.SCENARIO.FRESNEL(X) returns the
%   complex fresnel integral:
%
%              /X
%              |
%        Z  =  |  exp(1i * (pi/2)*s.^2) ds
%              |
%              /0

%#codegen

%   Copyright 2017-2021 The MathWorks, Inc.
%
%   Adapted from:
%      Cephes Math Library Release 2.8:  June, 2000
%      Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
%      Distributed under MIT License

if ~isreal(x)
    z = fresnelr(real(x)) + 1i*conj(fresnelr(-imag(x)));
else
    z = fresnelr(x);
end

function z = fresnelr(x)
% coder.noImplicitExpansionInFunction;
[sn, sd, cn, cd, fn, fd, gn, gd] = getConstants;

z = nan(size(x),'like',1i);

% compute positive argument
xabs = abs(x);
x2 = xabs .* xabs;

ismall = find(x2 < 2.5625);
x4 = x2(ismall);
x4 = x4 .* x4;
% Implicit expansion disabled here
z(ismall) = xabs(ismall) .* complex(  (((((cn(1).*x4 + cn(2)).*x4 + cn(3)).*x4 + cn(4)).*x4 + cn(5)).*x4 + cn(6)) ...
                                    ./ ((((((cd(1).*x4 + cd(2)).*x4 + cd(3)).*x4 + cd(4)).*x4 + cd(5)).*x4 + cd(6)).*x4 + cd(7)), ...
                      x2(ismall) .* (((((sn(1).*x4 + sn(2)).*x4 + sn(3)).*x4 + sn(4)).*x4 + sn(5)).*x4 + sn(6)) ...
                                    ./ ((((((sd(1).*x4 + sd(2)).*x4 + sd(3)).*x4 + sd(4)).*x4 + sd(5)).*x4 + sd(6)).*x4 + sd(7)));

% Asymptotic power series auxiliary functions for large argument 
ibig = find(2.5625 <= x2 & x2 <= 1367076676.0);
t = 1 ./ (pi * x2(ibig));
u = t .* t;
% Implicit expansion disabled here
dz = complex(u .* (((((((((  fn(1).*u + fn(2)).*u + fn(3)).*u + fn(4)).*u + fn(5)).*u + fn(6)).*u + fn(7)).*u + fn(8)).*u + fn(9)).*u + fn(10)) ...
               ./ (((((((((( fd(1).*u + fd(2)).*u + fd(3)).*u + fd(4)).*u + fd(5)).*u + fd(6)).*u + fd(7)).*u + fd(8)).*u + fd(9)).*u + fd(10)).*u + fd(11)), ...
             t .* (((((((((( gn(1).*u + gn(2)).*u + gn(3)).*u + gn(4)).*u + gn(5)).*u + gn(6)).*u + gn(7)).*u + gn(8)).*u + gn(9)).*u + gn(10)).*u + gn(11)) ...
               ./ (((((((((((gd(1).*u + gd(2)).*u + gd(3)).*u + gd(4)).*u + gd(5)).*u + gd(6)).*u + gd(7)).*u + gd(8)).*u + gd(9)).*u + gd(10)).*u + gd(11)).*u + gd(12))   ) - 1;

% Implicit expansion disabled here
z(ibig) = 0.5 + 0.5i + 1i*dz.*exp(0.5i*pi .* x2(ibig))./(pi .* xabs(ibig));

% clobber anything bigger
z(x2 > 1367076676.0) = 0.5 + 0.5i;

% odd symmetric for negative argument
z(x < 0.0) = -z(x < 0.0);

function [sn, sd, cn, cd, fn, fd, gn, gd] = getConstants
% The integrals are evaluated by a power series for x < 1.
% S(x) for small x 
sn = [ -2.99181919401019853726E3
        7.08840045257738576863E5
       -6.29741486205862506537E7
        2.54890880573376359104E9
       -4.42979518059697779103E10
        3.18016297876567817986E11];
sd = [  1.00000000000000000000E0
        2.81376268889994315696E2
        4.55847810806532581675E4
        5.17343888770096400730E6
        4.19320245898111231129E8
        2.24411795645340920940E10
        6.07366389490084639049E11];

% C(x) for small x 
cn = [ -4.98843114573573548651E-8
        9.50428062829859605134E-6
       -6.45191435683965050962E-4
        1.88843319396703850064E-2
       -2.05525900955013891793E-1
        9.99999999999999998822E-1];

cd = [  3.99982968972495980367E-12
        9.15439215774657478799E-10
        1.25001862479598821474E-7
        1.22262789024179030997E-5
        8.68029542941784300606E-4
        4.12142090722199792936E-2
        1.00000000000000000118E0];

% For x >= 1 auxiliary functions f(x) and g(x) are employed
% such that
%
% C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
% S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
%
% ACCURACY: Relative error.
%
% Arithmetic  function   domain     # trials      peak         rms
%   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
%   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
%

% Auxiliary function f(x) 
fn = [  4.21543555043677546506E-1
        1.43407919780758885261E-1
        1.15220955073585758835E-2
        3.45017939782574027900E-4
        4.63613749287867322088E-6
        3.05568983790257605827E-8
        1.02304514164907233465E-10
        1.72010743268161828879E-13
        1.34283276233062758925E-16
        3.76329711269987889006E-20];

fd = [  1.00000000000000000000E0
        7.51586398353378947175E-1
        1.16888925859191382142E-1
        6.44051526508858611005E-3
        1.55934409164153020873E-4
        1.84627567348930545870E-6
        1.12699224763999035261E-8
        3.60140029589371370404E-11
        5.88754533621578410010E-14
        4.52001434074129701496E-17
        1.25443237090011264384E-20];

% Auxiliary function g(x)
gn = [  5.04442073643383265887E-1
        1.97102833525523411709E-1
        1.87648584092575249293E-2
        6.84079380915393090172E-4
        1.15138826111884280931E-5
        9.82852443688422223854E-8
        4.45344415861750144738E-10
        1.08268041139020870318E-12
        1.37555460633261799868E-15
        8.36354435630677421531E-19
        1.86958710162783235106E-22];
        
gd = [  1.00000000000000000000E0
        1.47495759925128324529E0
        3.37748989120019970451E-1
        2.53603741420338795122E-2
        8.14679107184306179049E-4
        1.27545075667729118702E-5
        1.04314589657571990585E-7
        4.60680728146520428211E-10
        1.10273215066240270757E-12
        1.38796531259578871258E-15
        8.39158816283118707363E-19
        1.86958710162783236342E-22];
