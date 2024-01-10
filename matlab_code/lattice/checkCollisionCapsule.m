function [collisionFound, distance] = checkCollisionCapsule(p1, v1, D1, R1, p2, v2, D2, R2, exhaustive)
% This function is for internal use only. It may be removed in the future.
%checkCollisionCapsule Calculates the distance between line segments
%
%   COLLISIONFOUND = checkCollisionCapsule(P1,V1,D1,R1,P2,V2,D2,R2) performs
%   a pair-wise comparison of capsules in Set1 against capsules in Set2. If
%   the number of capsules in Set2 (e.g M) is an integer multiple of the number of
%   capsules in Set1 (e.g. N), then Set1 is compared against M/N horizontally
%   concatenated set of capsules in Set2. COLLISIONFOUND is returned as an 
%   N-by-(M/N) matrix of logicals, where the [i,j]-th element represents 
%   the comparison of the ith capsule in Set1 against the paired capsule 
%   in the jth concatenated set of capsules in Set2.
%
%       [A1 A2 .. AN] |^| [[B11 B12 .. B1N] [B21 B22 .. B2N] ... [...]]
% 
%                           -> [A1_B11 A1_B21     A1_BM1
%                               A2_B12 A2_B22     A2_BM2
%                                   ...
%                               AN_B1N AN_B2N ... AN_BMN] 
% 
%       An individual capsule is defined by a set of parameters:
%           P   - An xy or xyz point corresponding to the start of the capsule's axis
%           V   - An xy or xyz vector describing the direction of the capsule's axis
%           D   - The length of the line segment at the core of the capsule
%           R   - The capsule's radius
%
%                           ,.-+-------------------+-.,
%                          /       V                   \
%                Z  Y     (    P===>- - - - - - - -+    )
%                | /       \                        \R /
%                |/_ _X     `*-+-------------------+-*`
%                              |---------D---------|
%
%       For a set of capsules, P and V are DIM-by-NumCapsule matrices, where 
%       DIM is either 2 or 3 depending on whether the capsule is 2D (xy) or
%       3D (xyz). D and R can either be scalar, or NumCapsule-element vectors.
%
%   COLLISIONFOUND = checkCollisionCapsule(P1,V1,D1,R1,P2,V2,D2,R2,EXHAUSTIVE)
%   takes an optional logical flag, EXHAUSTIVE. If this is set to true,
%   then each capsule in Set1 is compared against ALL capsules in Set2. The
%   size of COLLISIONFOUND will be N-by-M, where the [i,j]-th element
%   represents the comparison of SetA's ith capsule with SetB's jth
%   capsule.
%
%           [A1 A2 .. AN] |^| [B1 B2 .. BN]
% 
%                           -> [A1_B1 A1_B2     A1_BM
%                               A2_B1 A2_B2     A2_BM
%                                   ...
%                               AN_B1 AN_B2 ... AN_BM]
%
%   [COLLISIONFOUND, DISTANCE] = checkCollisionCapsule(P1,V1,D1,R1,P2,V2,D2,R2,___)
%   returns the distance between pairs of capsules in Set1 and Set2. If a
%   pair of capsules are not in collision, then DISTANCE will contain the
%   distance between the nearest points on the capsules' boundary. If two
%   capsules are intersecting, the corresponding element of DISTANCE will 
%   contain NaN.
%
%   References:
%
%       [1] Distance calculation based on http://geomalgorithms.com/a07-_distance.html

%   Copyright 2020 The MathWorks, Inc.

%#codegen
    
    narginchk(8,9)
    
    numSeg1 = size(p1,2);
    numSeg2 = size(p2,2);
    dim = size(p1,1);
    
    if numSeg1 == 0
        collisionFound = false(0,1);
        distance = [];
        return;
    end
    
    if numSeg2 == 0
        collisionFound = false(numSeg1,1);
        if nargout == 2
            distance = inf(numSeg1,1);
        else
            distance = [];
        end
        return;
    end
    
    if nargin > 8 && exhaustive == true
        % The first set of line-segments will be checked exhaustively
        % against the second set
        
        % Calculate values needed to vectorize operations
        numChecks = numSeg1*numSeg2;
        
        % Replicate vectors
        v = repmat(v1.*repelem(D1,dim,1),1,numSeg2);
        u = repelem(v2.*repelem(D2,dim,1),1,numSeg1);
        
        % Calculate distance between beginning of line-segment pairs
        w0 = repelem(p2,1,numSeg1)-repmat(p1,1,numSeg2); % Vector from lineSeg1-base to lineSeg2-base
        
        % Calculate distance thresholds for each pair of radii
        combinedRSquared = (repmat(R1(:)',1,numChecks/numel(R1)) + repelem(R2(:)',1,numChecks/numel(R2))).^2;
        
        % Calculate height dimension of output
        outputHeight = numSeg1;
    else
        % Each line-segment in xy1's ith column will be checked against 
        % the set of line segments in xy2, corresponding to columns
        % i:size(xy1,2):end
        
        % Calculate the number of times set1 must be checked against set2
        numObj = numSeg2/numSeg1;
        validateattributes(numObj,{'numeric'},{'integer','positive'},'collisionCheckCapsuleVectorized','LineSegRatio');
        outputHeight = size(p1,2);

        % Replicate vectors
        v = repmat(v1.*D1,1,numObj);
        u = v2.*repelem(D2,dim,1);

        % Calculate distance between beginning of line-segment pairs
        w0 = p2-repmat(p1,1,numObj);

        % Calculate distance thresholds for each pair of radii
        combinedRSquared = (repmat(R1(:)',1,numObj*numSeg1/numel(R1))+R2).^2;

        numChecks = size(p1,2)*numObj;
    end
    
    a = sum(u.*u); % Unit Vector
    b = sum(v.*u);
    c = sum(v.*v); % Unit Vector
    d = sum(u.*w0);
    e = sum(v.*w0);
    
    % Find parametric values for nearest points on each line segment
    [tC, sC] = findNearestPoints(numChecks, a, b, c, d, e);
    
    % Calculate closest distance between pairs of line segments
    squareDist = sum((w0 + repmat(sC,dim,1).*u - repmat(tC,dim,1).*v).^2);
    
    % If distance is below the bounding radii of the two capsules, a
    % collision has occurred
    collisionFound = reshape(squareDist(:) <= combinedRSquared(:), outputHeight, []);
    
    if nargout == 2
        distance = reshape(sqrt(max(squareDist(:),0))-sqrt(combinedRSquared(:)), outputHeight, []);
        distance(distance<0) = nan;
    end
end

function [tC, sC] = findNearestPoints(numChecks, a,b,c,d,e)
%findNearestPoints Constrains the parametric variables to the length
%of each line segment to find the nearest points on each line

    % Calculate numerator and denominator for the parameters of the line-segments
    denom = a.*c - b.*b;
    tC = a.*e-b.*d;
    sC = b.*e-c.*d;
    
    % Precalc our floating-point check rather than repeatedly calculate
    sqrtEps = sqrt(eps);
    
    for i = 1:numChecks
        % Constrain parameter S of line 2 so that it falls within the
        % bounds of the line segment and recalculate T
        sD = denom(i);
        if sD < sqrtEps
            % Degenerate case where lines are collinear or parallel
            sC(i) = 0;
            sD = 1;
            tC(i) = e(i);
            tD = c(i);
        elseif sC(i) < 0
            sC(i) = 0;
            tC(i) = e(i);
            tD = c(i);
        elseif sC(i) > sD
            sC(i) = sD;
            tC(i) = e(i)+b(i);
            tD = c(i);
        else
            tD = denom(i);
        end
        
        % Constrain the parameter T such that it falls within the bounds of
        % line-seg1 while minimizing the distance between line 1 and line 2
        if tC(i) < 0
            tC(i) = 0;
            if -d(i) < 0
                sC(i) = 0;
            elseif -d(i) > a(i)
                sC(i) = sD;
            else
                sC(i) = -d(i);
                sD = a(i);
            end
        elseif tC(i) >= tD
            tC(i) = tD;
            if (-d(i) + b(i)) < 0
                sC(i) = 0;
            elseif (-d(i) + b(i) > a(i))
                sC(i) = sD;
            else
                sC(i) = -d(i) + b(i);
                sD = a(i);
            end
        end
        
        % Recalculate the parameters using the updated numerator and denominator
        if sC(i) < sqrtEps
            sC(i) = 0;
        else
            sC(i) = sC(i)/sD;
        end
        if tC(i) < sqrtEps
            tC(i) = 0;
        else
            tC(i) = tC(i)/tD;
        end
    end
end