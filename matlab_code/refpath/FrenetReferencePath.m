classdef FrenetReferencePath < handle
%This class is for internal use only. It may be removed in the future.

%FrenetReferencePath Generate a reference path using waypoints
%   FrenetReferencePath It fits a parametric curve over the waypoints
%   and generate a discrete set of states which contains x, y, theta,
%   kappa, dkappa and s
    properties(Access=public)
        PathPoints
        Waypoints
        MaxNumWaypoints
        %Length Total length of the path
        Length
        %PathManager Class responsible for managing waypoint data
        PathManager
        DiscretizationDistance
        BoundingBoxes
    end
    methods
        function obj = FrenetReferencePath(waypoints, varargin)
            nvPairs.DiscretizationDistance = 0.05;
            nvPairs.MaxNumWaypoints = 10000;
            % Set discretization step size
            obj.DiscretizationDistance = nvPairs.DiscretizationDistance;
            obj.Waypoints = waypoints;
            obj.generateStates(waypoints);
        end

        function pathPoints = interpolate(obj, arcLength)
            pathPoints = pathInterpolate(obj.PathManager.SegStarts, arcLength(:));
        end
        
        function [pts, dists, segIdx, ptIdx] = nearestPoint(obj, xyPathPoints)
             [pts,dists,segIdx,ptIdx] = nearestPointIterativeImp(xyPathPoints(:,1:2), obj.PathManager.SegStarts, obj.BoundingBoxes, obj.Length);
        end

        function pathPoints = closestPoint(obj, xyPointss)
            %closestPoint Find the closest point along path from (x,y) coordinate
            %   Find the closest path point on the reference path from the
            %   input point (x,y)
            ss = 0:0.1:obj.Length;
            obj.PathPoints = obj.interpolate(ss);
            pathPoints = zeros(size(xyPointss,1), 6);
            for i = 1:size(xyPointss,1)
                xyPoints = xyPointss(i,:);
                dist = sqrt((xyPoints(1) - obj.PathPoints(:,1)).^2 +(xyPoints(2) - obj.PathPoints(:,2)).^2);
                [~,minIdx] = min(dist);
                if minIdx == 1
                    cur =  minIdx;
                    next = minIdx + 1;
                    dxy1 =  obj.PathPoints(cur,1:2) - xyPoints;
                    dxy2 =  obj.PathPoints(next,1:2) - obj.PathPoints(cur,1:2);
                    dxy2 = dxy2/norm(dxy2);
                    dotval = dxy1*(dxy2');
                    if dotval > 0
                        pathPoint = obj.PathPoints(cur,:);
                    else
                        pathPoint = obj.interpolate(-dotval);
                    end
                elseif minIdx == length(ss)
                    cur = minIdx - 1;
                    next = minIdx;
                    dxy1 =  xyPoints - obj.PathPoints(next,1:2);
                    dxy2 =  obj.PathPoints(next,1:2) - obj.PathPoints(cur,1:2);
                    dxy2 = dxy2/norm(dxy2);
                    dotval = dxy1*(dxy2');
                    if dotval > 0
                        pathPoint = obj.PathPoints(next,:);
                    else
                        pathPoint = obj.interpolate(obj.PathPoints(end,6)+dotval);
                    end
                else
                    % TODO 3点投影关系
                    pathPoint = obj.PathPoints(minIdx,:);
                end
                pathPoints(i, :) = pathPoint;
            end
        end
        
        function frenetState = global2frenet(obj, globalStates)
        %global2frenet Convert global states to Frenet states
        % Get path point on reference path closest to given input state
%             refPathPoints = obj.closestPoint(globalStates(:,1:2));
            refPathPoints = obj.nearestPoint(globalStates(:,1:2));
            % Transform to Frenet states using path point as reference
            frenetState = CartesianFrenetConversions.cartesian2Frenet(refPathPoints, globalStates);
        end
        
        function globalStates = frenet2global(obj, frenetStates)
        %frenet2global Convert Frenet states to global states
        % Get path point on reference path closest to given input state
            refPathPoint = obj.interpolate(frenetStates(:,1));
            globalStates = CartesianFrenetConversions.frenet2Cartesian(refPathPoint, frenetStates);
        end

    end
    
    methods (Access = public)
        function generateStates(obj, waypoints)
            %generateStates Generate path points from waypoints
            %   generateStates fits a clothoid on the
            %   waypoints and computes poses, curvature and curvature
            %   derivative
            % Fit a clothoid over waypoints
            switch size(waypoints,2)
                case 2
                    xy = waypoints;
                    obj.fitClothoid(xy);
                case 3
                    xy = waypoints(:,1:2);
                    th = waypoints(:,3);
                    obj.fitClothoid(xy, th);
                otherwise
                    coder.internal.error('nav:navalgs:referencepathfrenet:IncorrectNumCols')
            end
            obj.BoundingBoxes = obj.calculateBoundingBoxes(waypoints);
        end
        
        function bounds = calculateBoundingBoxes(obj, wp)
        %calculateBoundingBoxes Calculates the AABB limits for each section of the curve
            numSeg = size(wp,1)-1;
            bounds = zeros(4,numSeg);
            ss = obj.PathManager.SegStarts;
            for i = 1:numSeg
                bounds(:,i) = findClothoidBoundingBox(ss(i,1),ss(i,2),ss(i,3),ss(i,4),ss(i,5),obj.PathManager.Arclengths(i));
            end
        end
        
        function fitClothoid(obj, waypoints, courseInput)
        %fitClothoid Fits a clothoid on the waypoints and computes states
            if nargin == 2
                % clothoidG2fitCourse calc waypoints headangle
                course = clothoidG2fitCourse(waypoints);
            else
                course = courseInput;
            end
            n = size(waypoints,1);
            % Obtain the initial positions
            initialPosition = complex(waypoints(:,1), waypoints(:,2));
            arcLengths = zeros(n,1);
            % Calculate starting curvature, final curvature, and length of each segment.
            [initialCurvature, finalCurvature, arcLengths(1:end-1)] = clothoidG1fit2(initialPosition(1:n-1),course(1:n-1),initialPosition(2:n),course(2:n));
            % Store the path conditions at each waypoint and the arclength
            % of each clothoid segment.
            cumulativeLength = cumsum(circshift(arcLengths,1));
            segStarts = [real(initialPosition) imag(initialPosition) robotics.internal.wrapToPi(course) [initialCurvature; finalCurvature(end)] [(finalCurvature-initialCurvature)./arcLengths(1:end-1); (finalCurvature(end)-initialCurvature(end))/arcLengths(end-1)] cumulativeLength];
            obj.PathManager.Arclengths     = arcLengths;
            obj.PathManager.SegStarts      = segStarts;
            obj.Length         = cumulativeLength(end);
        end
    end
   
    methods (Hidden, Static)
        function [minDist, minIdx, segDist] = nearestBoundingBox(bounds, xy)
        %nearestBoundingBox Find nearest AABB to point
            numSeg = size(bounds,2);
            segDist = inf(1,numSeg);
            
            % Find distance to each boundingbox
            for i = 1:numSeg
                if xy(1) < bounds(1,i)
                    % Left-side
                    xDist = xy(1)-bounds(1,i);
                elseif xy(1) > bounds(2,i)
                    % Right-side
                    xDist = xy(1)-bounds(2,i);
                else
                    % Between horizontal
                    xDist = 0;
                end
                if xy(2) < bounds(3,i)
                    % Bottom-side
                    yDist = xy(2)-bounds(3,i);
                elseif xy(2) > bounds(4,i)
                    % Top-side
                    yDist = xy(2)-bounds(4,i);
                else
                    % Between vertical
                    yDist = 0;
                end
                
                % Store distance to AABB
                segDist(i) = (xDist^2+yDist^2)^(1/2);
            end
            
            % Return nearest BB index
            [minDist,minIdx] = min(segDist);
        end
        
        function [x1,y1,th1,k1] = clothoid(x0,y0,th0,k0,dk,L)
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
    end
end
