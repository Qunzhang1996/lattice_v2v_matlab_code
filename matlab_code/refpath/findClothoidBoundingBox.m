function bounds = findClothoidBoundingBox(x0,y0,th0,k0,dk,L)
%This function is for internal use only. It may be removed in the future.

%findClothoidBoundingBox Finds the axially aligned bounds for the provided clothoid segment
%
%   Finds the xy limits that bound a segment along a curve. The beginning of
%   the segment is located at (x0,y0), with tangent-angle and curvature of
%   th0 and k0, respectively. The change in curvature, dk, is constant, and
%   the arclength of the segment is L. All inputs must be scalar, and L
%   must be non-negative.
%
%   Note:
%       - If k0,dk are BOTH zero, the curve represents a line-segment.
%       - If dk is zero, the curve represents an arc.
%       - If dk is non-zero, the curve represents a clothoid.
%
%   Copyright 2020 The MathWorks, Inc.

%#codegen

    validateattributes(L,{'numeric'},{'nonnegative','vector'},'findClothoidBoundingBox','L',6)
    if dk == 0
        if k0 == 0
            bounds = calculateBoundingBoxLine(x0,y0,th0,L);
        else
            bounds = calculateBoundingBoxArc(x0,y0,th0,k0,L);
        end
    else
        bounds = calculateBoundingBoxClothoid(x0,y0,th0,k0,dk,L);
    end
end

function bounds = calculateBoundingBoxLine(x0,y0,th0,L)
%calculateBoundingBoxLine Finds the axially-aligned limits bounding a line-segment

    coder.inline('always')
    
    x = x0+cos(th0)*L;
    y = y0+sin(th0)*L;
    bounds = [min(x,x0);max(x,x0);min(y,y0);max(y,y0)];
end

function bounds = calculateBoundingBoxArc(x0,y0,th0,k0,L)
%calculateBoundingBoxArc Finds the axially-aligned limits bounding an arc

    coder.inline('always');

    % Calculate properties of the arc
    dTh = k0*L;
    th1 = th0+dTh;
    R = 1/k0;
    xCenter = x0-R*sin(th0);
    yCenter = y0+R*cos(th0);

    if dTh >= 2*pi
    %Full circle
        r = abs(R);
        bounds = [xCenter-r; xCenter+r; yCenter-r; yCenter+r];
    else
    %Arc segment
        % Calculate the tangent angle of the first critical point 
        % encountered when following the segment from s0->s0+L.
        if k0 > 0
        %Positive curvature -> counter-clockwise motion
            thInit = ceil(th0/(pi/2))*pi/2;
        else
        %Negative curvature -> clockwise motion
            thInit = floor(th0/(pi/2))*pi/2;
        end
        
        % If the endpoint is encountered before the first critical point,
        % both points lie in the same quadrant of the circle.
        if (dTh-(thInit-th0))*sign(k0) <= 0
            %Segment tangent-angles lie in the same quadrant
            th = th1;
        else
            %Segment tangent-angles lie in separate quadrants, critical points
            %must be considered.
            th = [thInit:(pi/2)*sign(k0):th1 th1];
        end
        xVals = [x0 xCenter+R*sin(th)];
        yVals = [y0 yCenter-R*cos(th)];
        bounds = [min(xVals); max(xVals); min(yVals); max(yVals)];
    end
end

function bounds = calculateBoundingBoxClothoid(x0,y0,th0,k0,dk,L)
%calculateBoundingBoxClothoid Finds the axially-aligned limits bounding a clothoid
%
%   Given the input arguments that describe a clothoid, this subroutine 
%   attempts to find all critical points (arclengths at which the curve's 
%   tangent-angle -> [0,pi,pi/2,-pi/2]*i, i=..,-1,0,1,..) lying between the
%   segment boundaries.
%   
%   To simplify the problem, the subroutine:
%
%   1) Finds the arclength between (x0,y0) and the clothoid's inflection point (s0)
%   2) Determines whether the segment:
%       a) Excludes the inflection point -> Immediately solvable
%       b) Contains the inflection point
%           - Breaks the segment into two pieces that lie on opposite sides 
%             of the inflection
%           - Flip the inward-facing segment so that they both point away
%           - Finds AABB for each piece
%           - Returns the union of both AABB

    % Calculate arclength of base-point relative to inflection
    s0 = k0/dk;
    sSign = k0*dk;
    
    if sSign >= 0
    % If k0,dk are of the same sign, then segment points away from inflection pt
        bounds = findLargestBoundingBox(x0,y0,th0,k0,dk,s0,s0+L);
    else
    % Segment points towards inflection
        if L*dk*dk + sSign <= eps
        % Segment does not cross inflection point
            [xp,yp,thp,kp] = FrenetReferencePath.clothoid(x0,y0,th0,k0,dk,L);
            bounds = findLargestBoundingBox(xp,yp,thp+pi,-kp,dk,-s0-L,-s0);
        else
        % Inflection point contained between [s0, s0+L]
        % Search both sections and return the bounding box of both

            % Find center point
            [xBar,yBar,thBar] = FrenetReferencePath.clothoid(x0,y0,th0,k0,dk,-s0);

            % Find bounds for each segment
            bounds1 = calculateBoundingBoxClothoid(xBar,yBar,thBar+pi,0,dk,abs(s0));
            bounds2 = calculateBoundingBoxClothoid(xBar,yBar,thBar,0,dk,L+s0);

            % Merge boxes
            bounds = zeros(4,1);
            bounds([1;3;2;4]) = [min(bounds1([1 3]),bounds2([1 3]));max(bounds1([2 4]),bounds2([2 4]))];
        end
    end
end

function bounds = findLargestBoundingBox(x0,y0,th0,k0,dk,s0,sMax)
%findLargestBoundingBox Finds the AABB that bounds all critical/endpoints on the interval [s0 sMax]

    % "Standardize" segment such that we only search for positive roots
    thNew = th0*sign(dk);
    kNew = k0*sign(dk);
    DK = abs(dk);
    thBar = DK/2*s0^2-kNew*s0+thNew;
    
    % Calculate all feasible critical points
    thetaIdxAll = fix((thBar+s0^2*DK/2)/(pi/2))+[0 1 2 3 4];
    rtAll = -2*DK*(thBar-thetaIdxAll*(pi/2));
    sCAll = sqrt(rtAll(rtAll>0))/DK;
    
    % Exclude points that lie outside the interval [s0 sMax]
    s = [sCAll(sCAll > s0 & sCAll < sMax) sMax]-s0;
    
    [x,y] = FrenetReferencePath.clothoid(x0,y0,th0,k0,dk,s);
    
    xVals = [x0 x];
    yVals = [y0 y];
    
    bounds = [min(xVals);max(xVals);min(yVals);max(yVals)];
end