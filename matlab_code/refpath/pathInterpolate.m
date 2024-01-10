function pathPoints = pathInterpolate(segStarts, arcLength)
% This class is for internal use only. It may be removed in the future.

%interpolate Evaluate the piecewise continous clothoid at provided arclengths
%
%   Takes in an S-segment piecewise-continuous (C2) clothoid-spline,
%   where each segment is defined by a row in SEGSTARTS, taking the form 
%   [x0,y0,th0,k0,dk,S]. (S is relative to the beginning of the PWC curve). 
%   ARCLENGTH is an N-element vector of arclengths at which the
%   clothoid-spline should be evaluated.
%
%   Example:
%       % Create reference path and extract state of the path at each waypoint.
%       waypoints = [0 0 0; 10 5 pi/2; -10 0 -pi/2];
%       refPath = referencePathFrenet([0 0 0; 10 5 pi/2; -10 0 -pi/2]);
%       segStarts = refPath.closestPoint(waypoints(:,1:2));
%
%       % Interpolate along the path at evenly spaced intervals
%       lengths = linspace(0,segStarts(end),100);
%       interpolatedPts = nav.algs.internal.impl.pathInterpolate(segStarts,lengths);
%
%       % Display discretized curve
%       plot(interpolatedPts(:,1),interpolatedPts(:,1));

%   Copyright 2020 The MathWorks, Inc.

    %#codegen
    
    % Find the clothoid to which each arclength belongs
    segIdx = discretize(arcLength, [-inf; segStarts(2:end,end); inf]);
    
    % Allocate output
    endPt = segStarts(end,:);
    pathPoints = repmat(endPt,numel(arcLength),1);
    S = max(arcLength,0);
    for i = 1:numel(segIdx)
        idx = segIdx(i);
        if S(i) < endPt(end)
            if S(i) == 0
                pathPoints(i,:) = segStarts(1,:);
            else
                [x,y,th,k] = FrenetReferencePath.clothoid(segStarts(idx,1),segStarts(idx,2),segStarts(idx,3),segStarts(idx,4),segStarts(idx,5),S(i)-segStarts(idx,6));
                pathPoints(i,:) = [x y robotics.internal.wrapToPi(th) k segStarts(idx,5) S(i)];
            end
        end
    end
end