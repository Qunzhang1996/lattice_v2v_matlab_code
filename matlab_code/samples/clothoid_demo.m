%clothoid_demo.m

L = 0:0.01:25;
[xs, ys,~,ks] = clothoidF(0, 0, 0, 0, 0.1, L);
plot(xs, ys, 'marker', '.')
axis equal;

figure
plot(L,ks, 'marker', '.')
title('曲率')

function [x1,y1,th1,k1] = clothoidF(x0,y0,th0,k0,dk,L)
%clothoid Construct terminal points for clothoid
if isscalar(L) && L(1) == 0
    x1 = x0; y1 = y0; th1 = th0; k1 = k0;
else
    k1  = k0+dk*L;
    th1 = dk/2*L.^2 + k0*L + th0;
    % Calc integrated xy
    xy = fresnelg(L, dk, k0, th0);
    % Add to base location
    x1 = x0 + real(xy);
    y1 = y0 + imag(xy);
end
end

