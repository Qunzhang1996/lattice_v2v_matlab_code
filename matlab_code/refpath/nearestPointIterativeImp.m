function [pts, dist, segIdx, ptIdx] = nearestPointIterativeImp(xyPts, segStarts, bounds, sMax)
% This class is for internal use only. It may be removed in the future.

%nearestPointIterative Returns closest point along a piecewise continuous clothoid
%
%   This function is given an Nx2 set of query points, XYPTS, and the
%   definition of an M-segment piecewise continuous clothoid, SEGSTARTS. It
%   is also provided with BOUNDS, a 4xM matrix of XY limits, which define
%   axially-aligned boxes which bound the curve in each segment. Lastly, it
%   is passed the maximum arclength, sMax, which defines the length of the
%   final segment.
%
%   [PTS, DISTS, SEGIDX, PTIDX] = nearestPtIterative(XYPTS, PATHPTS, BOUNDS, SEGMENTINDICES, NUMSEG)
%       Inputs:
%           XYPTS     - Set of XY points (Nx2)
%           SEGSTARTS - M-by-[x y O k dk S] definition of clothoid at start of each PWC clothoid segment
%           BOUNDS    - [Left; Right; Bottom; Top]-by-M limits for boxes that bound each PWC clothoid segment
%           SMAX      - Total length of path        
% 
%       Outputs:
%           PTS    - Nearest point on the path (N-by-[x y theta kappa dkappa s])
%           DIST   - Distance between point and path (Nx1)
%           SEGIDX - Clothoid segment index containing PTS (Nx1)
%           PTIDX  - Index of nearest discretized point along entire path
%
%   Finding the closest point is broken into four operations. For each query point:
%
%       1) Find the closest bounding box (nearestBoundingBox)
%       2) Find the true min distance to bounded segment (nearestPointOnClothoid)
%       3) Loop through all remaining bounding box. If their bounds fall 
%          within the current min dist, find the true min dist of that
%          segment and retain the smallest.
%       4) Calculate/return the distance and point along the minimum clothoid 
%          using the sub-segment arclength, sMin.

%   Copyright 2020-2021 The MathWorks, Inc.

    %#codegen
    numPts = size(xyPts,1);
    pts    = zeros(numPts,6);
    dist  = zeros(numPts,6);
    segIdx = zeros(numPts,1);
    ptIdx  = zeros(numPts,1);
    
    for i = 1:numPts
        xy = xyPts(i,:);
        % Find distance to all AABB and return the nearest AABB
        [~, minBBIdx, allDist] = FrenetReferencePath.nearestBoundingBox(bounds, xy);

        cInfoBest = segStarts(minBBIdx,:);
        if minBBIdx == size(segStarts,1)
            dS = sMax-cInfoBest(end);
        else
            dS = segStarts(minBBIdx+1,6)-cInfoBest(end);
        end

        [sMin, minDist] = nearestPointOnClothoid(xy,cInfoBest(1),cInfoBest(2),cInfoBest(3),cInfoBest(4),cInfoBest(5),dS);
        sIdx = minBBIdx;
        
        for j = 1:size(segStarts,1)-1
            if j ~= minBBIdx && minDist >= allDist(j)
                % Distance between query pt and BB is less than
                % closest point on current segment
                
                % Analytic method
                cInfo = segStarts(j,:);
                if j == size(segStarts,1)
                    dS = sMax-cInfo(end);
                else
                    dS = segStarts(j+1,6)-cInfo(end);
                end
                [s, d] = nearestPointOnClothoid(xy,cInfo(1),cInfo(2),cInfo(3),cInfo(4),cInfo(5),dS);

                if d < minDist
                    cInfoBest(1,:) = reshape(cInfo,1,6);
                    minDist = d;
                    sMin = s;
                    sIdx = j;
                end
            end
        end

        % Add to output results
        [x,y,th,k] = FrenetReferencePath.clothoid(cInfoBest(1,1),cInfoBest(1,2),cInfoBest(1,3),cInfoBest(1,4),cInfoBest(1,5), sMin);
        S          = cInfoBest(end)+sMin;
        pts(i,:)   = [x y th k cInfoBest(5) S];
        dist(i)    = minDist;
        segIdx(i)  = sIdx;
    end
end