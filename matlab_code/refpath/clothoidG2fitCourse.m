function course = clothoidG2fitCourse(waypoints)
%MATLABSHARED.TRACKING.INTERNAL.SCENARIO.CLOTHOIDG2FITCOURSE find tangent angles at waypoints for G2
%clothoid fit approximate course with discrete clothoid fit upsampled 1024x
%
%   This function is for internal use only and may be removed in a later
%   release.
%

%   Copyright 2017-2019 The MathWorks, Inc.

%#codegen
n = size(waypoints,1);
course = zeros(n,1);
if n == 2
    % When there are only 2 waypoints, course is computed directly without calling dclothoid
    course(1) = angle(complex(waypoints(2,1)-waypoints(1,1),waypoints(2,2)-waypoints(1,2)));
    course(n) = course(1);
else
    % Get initial approximation to course angles
    up=1024;
    [u,v] = dclothoid(waypoints(:,1),waypoints(:,2));
    course(1) = angle(complex(u(2)-u(1),v(2)-v(1)));
    course(n) = angle(complex(u(end)-u(end-1),v(end)-v(end-1)));
    for i=1:n-2
        course(i+1) = angle(complex(u(i*up+1)-u(i*up),v(i*up+1)-v(i*up)));
    end
    courselsq = LMFsolve(@(x) clothresid(waypoints,x), course);
    course = courselsq;
end

function [kerr,Jtri] = clothresid(waypoints, course)
% obtain the (horizontal) initial positions
hip = complex(waypoints(:,1), waypoints(:,2));

% get the total number of control points
n = size(waypoints(:,1),1);

% for each segment, extract the initial and final curvatures and their
% partial derivatives with respect to the initial and final course angles
[k0, k1, ~, dk0_dc0, dk0_dc1, dk1_dc0, dk1_dc1] = clothoidG1fit2(hip(1:n-1),course(1:n-1),hip(2:n),course(2:n));

% Set desired curvature at endpoints
k0_begin = 0;
k1_end = 0;

% return the differences in curvature
kerr = [k0_begin-k0(1); k1(1:end-1)-k0(2:end); k1(end)-k1_end];

% The Jacobian is an N by N tridiagonal matrix
upper_diag = -dk0_dc1;
lower_diag =  dk1_dc0;

center_diag = [k0_begin       - dk0_dc0(1); 
               dk1_dc1(1:n-2) - dk0_dc0(2:n-1); 
               dk1_dc1(n-1)   - k1_end];
           
% J = sparse([1:n-1 1:n 2:n],[2:n 1:n 1:n-1],[upper_diag; center_diag; lower_diag],n,n);
Jtri = [[lower_diag; 0], center_diag, [upper_diag; 0]];