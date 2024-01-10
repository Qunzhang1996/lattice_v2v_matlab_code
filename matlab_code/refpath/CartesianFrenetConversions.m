classdef CartesianFrenetConversions
%This class is for internal use only. It may be removed in the future.

%CartesianFrenetConversions Conversion functions for the two frames
%
%   CartesianFrenetConversions properties:
%
%       cartesian2Frenet                    - Cartesian to Frenet
%                                             conversion method
%
%       frenet2Cartesian                    - Frenet to Cartesian
%                                             conversion method

%   Copyright 2019 The MathWorks, Inc.
%#codegen
    methods(Static)
        function frenetState = cartesian2Frenet(refState, queryState)
        %cartesian2Frenet convert from Cartesian to Frenet
        %   frenetState converts states [x, y, theta, kappa, v, a] to
        %   [s, sDot, sDotDot, l, lPrime, lPrimePrime] using the
        %   reference point which contains [x,y,theta,kappa,dkappa,s]

        % Validate number of input arguments
            narginchk(2,2)

            % Segregate reference states
            refS = refState(:,6);
            refX = refState(:,1);
            refY = refState(:,2);
            refTheta = robotics.internal.wrapToPi(refState(:,3));
            refKappa = refState(:,4);
            refKappaPrime = refState(:,5);

            % Segregate query Cartesian states
            x = queryState(:,1);
            y = queryState(:,2);
            v = queryState(:,5);
            a = queryState(:,6);
            theta = robotics.internal.wrapToPi(queryState(:,3));
            kappa = queryState(:,4);

            dx = x - refX;
            dy = y - refY;

            refCosTheta = cos(refTheta);
            refSinTheta = sin(refTheta);

            % Normal at the root point
            refNormal = refCosTheta .* dy - refSinTheta .* dx;
            % Compute lateral deviation
            l = vecnorm([dx, dy]')' .* sign(refNormal);

            deltaTheta = robotics.internal.angdiff(refTheta,theta);

            tanDeltaTheta = tan(deltaTheta);
            cosDeltaTheta = cos(deltaTheta);

            oneMinusKappaRefL = 1 - refKappa .* l;

            % Throw an error at extreme curvature or deltaTheta >= pi/2
            coder.internal.errorIf(isequal(or(oneMinusKappaRefL <= 0,...
                                              abs(deltaTheta) >= pi/2),ones(1,size(deltaTheta,2))), ...
                                   'shared_autonomous:cartesianFrenetConversions:singularity');

            % Compute derivative of lateral deviation w.r.t. arc length
            lPrime = oneMinusKappaRefL .* tanDeltaTheta;

            kappaRefLPrime = refKappaPrime .* l + refKappa .* lPrime;

            % Compute second derivative of lateral deviation w.r.t. arc length
            lPrimePrime = -kappaRefLPrime .* tanDeltaTheta + oneMinusKappaRefL ./ cosDeltaTheta ./ cosDeltaTheta .* (kappa .* oneMinusKappaRefL ./ cosDeltaTheta - refKappa);

            % Arc length is same as reference (root point) arc length
            s = refS;

            % Velocity in lateral direction
            sDot = v .* cosDeltaTheta ./ oneMinusKappaRefL;

            deltaThetaPrime = oneMinusKappaRefL ./ cosDeltaTheta .* kappa - refKappa;

            % Acceleration in lateral direction
            sDotDot = (a .* cosDeltaTheta - sDot.^2 .* (lPrime .* deltaThetaPrime - kappaRefLPrime)) ./ oneMinusKappaRefL;

            % Consolidate all the states as a 1x6 vector
            frenetState = [s sDot sDotDot l lPrime lPrimePrime];

            % Remove residuals by rounding off to order of sqrt(eps). See g2049087
            roundOrder = ceil(abs(log(sqrt(eps))/log(10)));
            frenetState = round(frenetState.*10^roundOrder)./10^roundOrder;
        end

        function cartesianStates = frenet2Cartesian(refState, queryState)
        %frenet2Cartesian Convert from Frenet to Cartesian
        %   frenetState converts [s, sDot, sDotDot, l, lPrime, lPrimePrime]
        %   to [x, y, theta, kappa, v, a] using the reference point which contains
        %   [x, y, theta, kappa, dkappa, s]

        % Validate number of input arguments
            narginchk(2,2)

            % Segregate reference states
            refX = refState(:,1);
            refY = refState(:,2);
            refTheta = robotics.internal.wrapToPi(refState(:,3));
            refKappa = refState(:,4);
            refKappaPrime = refState(:,5);

            % Segregate query Frenet states
            sDot = queryState(:,2);
            sDotDot = queryState(:,3);
            l = queryState(:,4);
            lPrime = queryState(:,5);
            lPrimePrime = queryState(:,6);

            refCosTheta = cos(refTheta);
            refSinTheta = sin(refTheta);

            % Compute x
            x = refX - refSinTheta .* l;

            % Compute y
            y = refY + refCosTheta .* l;

            oneMinusKappaRefL = 1 - refKappa .* l;

            tanDeltaTheta = lPrime ./ oneMinusKappaRefL;
            deltaTheta = atan2(lPrime, oneMinusKappaRefL);
            cosDeltaTheta = cos(deltaTheta);

            % Throw an error at extreme curvature or deltaTheta >= pi/2
            coder.internal.errorIf(isequal(or(oneMinusKappaRefL <= 0,...
                                              abs(deltaTheta) >= pi/2),ones(1,size(deltaTheta,2))), ...
                                   'shared_autonomous:cartesianFrenetConversions:singularity');

            % Compute theta
            theta = robotics.internal.wrapToPi(deltaTheta + refTheta + 2*pi);

            kappaRefLPrime = refKappaPrime .* l + refKappa .* lPrime;

            % Compute kappa
            kappa = (((lPrimePrime + kappaRefLPrime .* tanDeltaTheta) .* cosDeltaTheta.^2) ./ oneMinusKappaRefL + refKappa) .* cosDeltaTheta ./ oneMinusKappaRefL;

            % Compute speed in the direction of vehicle heading
            v = sDot .* (oneMinusKappaRefL ./ cosDeltaTheta);

            deltaThetaPrime = oneMinusKappaRefL ./ cosDeltaTheta .* kappa - refKappa;

            % Compute acceleration in the direction of vehicle heading
            a = (sDotDot .* oneMinusKappaRefL ./ cosDeltaTheta) + (sDot.^2 ./ cosDeltaTheta) .* (lPrime .* deltaThetaPrime - kappaRefLPrime);

            % Consolidate all the states as a 1x6 vector
            cartesianStates = [x y theta kappa v a];

            % Remove residuals by rounding off to order of sqrt(eps). See g2049087
            roundOrder = ceil(abs(log(sqrt(eps))/log(10)));
            cartesianStates = round(cartesianStates.*10^roundOrder)./10^roundOrder;
        end
    end
end
