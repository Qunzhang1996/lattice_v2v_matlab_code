classdef QuinticQuarticTrajectory
%This class is for internal use only. It may be removed in the future.

%QuinticQuarticTrajectory Generate 4th or 5th order polynomial
%
%   QuinticQuarticTrajectory methods:
%       QuinticQuarticTrajectory        - Constructor which takes in
%                                         variable value P and start and
%                                         end conditions at variable = 0
%                                         and variable = P
%
%       evaluate                        - Evaluate Nth order derivative
%                                         value of the polynomial at
%                                         the given variable value P
%

%   Copyright 2019 The MathWorks, Inc.

%#codegen
    properties
        %Coefficients - Coefficients the quintic polynomial
        %   The 6-vector of coefficients C used to define a quintic or a
        %   quartic polynomial Y in terms of variable X such that
        %   Y = C(6)*X^5 + C(5)*X^(4) + ... + C(2)*X + C(1)
        Coefficients

        %StartCondition - Initial conditions of the polynomial
        %   The 3-vector containing 0th, 1st, and 2nd order derivatives of
        %   the polynomial at variable X = 0
        StartCondition

        %EndCondition - Final conditions of the polynomial
        %   The 3-vector containing 0th, 1st, and 2nd order derivatives of
        %   the polynomial at variable X = Variable
        EndCondition

        %Variable - Indeterminate of the quintic function i.e Y(Variable)
        Variable
    end

    methods
        function obj = QuinticQuarticTrajectory(startCondition, endCondition, variable)
        %QuinticQuarticTrajectory Constructor for initializing the class object

        % Validate number of input arguments
            narginchk(3,3)

            obj.StartCondition = startCondition;
            obj.EndCondition = endCondition;
            obj.Variable = variable;
            if isnan(endCondition(1))
                obj.Coefficients = obj.computeQuarticTrajectory();
            else
                obj.Coefficients = obj.computeQuinticTrajectory();
            end
        end

        function result = evaluate(obj, order, variable)
        %evaluate Evaluate nth derivative of polynomial at variable

            validateattributes(order, {'double'},{'nonempty','integer','finite','<=',5,'>=',0,'vector'},'evaluate','order');

            result = zeros(length(order),length(variable));
            for i=1:length(order)
                currentOrder = order(i);
                switch currentOrder

                    % Evaluate quintic polynomial at variable
                  case 0
                    result(i,:) =  ((((obj.Coefficients(6) .* variable + obj.Coefficients(5)) .* variable + obj.Coefficients(4)) .* variable + obj.Coefficients(3)) .* variable + obj.Coefficients(2)) .* variable + obj.Coefficients(1);

                    % Evaluate first order derivative of quintic polynomial at variable
                  case 1
                    result(i,:) =  (((5.0 * obj.Coefficients(6) .* variable + 4.0 * obj.Coefficients(5)) .* variable + 3.0 * obj.Coefficients(4)) .* variable + 2.0 * obj.Coefficients(3)) .* variable + obj.Coefficients(2);

                    % Evaluate second order derivative of quintic polynomial at variable
                  case 2
                    result(i,:) =  (((20.0 * obj.Coefficients(6) .* variable + 12.0 * obj.Coefficients(5)) .* variable) + 6.0 * obj.Coefficients(4)) .* variable + 2.0 * obj.Coefficients(3);

                    % Evaluate third order derivative of quintic polynomial at variable
                  case 3
                    result(i,:) =  (60.0 * obj.Coefficients(6) .* variable + 24.0 * obj.Coefficients(5)) .* variable + 6.0 * obj.Coefficients(4);

                    % Evaluate forth order derivative of quintic polynomial at variable
                  case 4
                    result(i,:) =  120.0 * obj.Coefficients(6) .* variable + 24.0 * obj.Coefficients(5);

                    % Evaluate fifth order derivative of quintic polynomial at variable
                  case 5
                    result(i,:) =  120.0 * obj.Coefficients(6);
                end
            end

        end
    end

    methods(Access=private)
        function coefficients = computeQuinticTrajectory(obj)
        %computeTrajectory Computes the coefficients of the polynomial
        %   This function computes the quintic polynomial coefficients
        %   given the zeroth, first, and second derivative at P = 0 and
        %   P = finalVariable

        % Code taken from quinticpolytraj generateQuinticCoeffs method

        % First three coefficients can be found by setting P=0 using
        % boundary conditions
            coeffVec = [obj.StartCondition(1) obj.StartCondition(2) obj.StartCondition(3)/2.0 0 0 0]';

            % The last three coefficients may be found by substituting the first three
            % into the expressions for zeroth, first, and second derivative at
            % P=finalVariable:
            % [f(P) fd(P) fdd(P)]' = matP0*[C1 C2 C3]' + matPFinal*[C4 C5 C6]
            % where matPFinal =
            %   [finalVariable^3    finalVariable^4     finalVariable^5; ...
            %    3*finalVariable^2  4*finalVariable^3   5*finalVariable^4; ...
            %    6*finalVariable    12*finalVariable^2  20*finalVariable^3];
            % This expression may be solve for C4, C5, and C6

            matP0 = [1 obj.Variable obj.Variable^2; 0 1 2*obj.Variable; 0 0 2];

            B = [obj.EndCondition(1); obj.EndCondition(2); obj.EndCondition(3)] - matP0*coeffVec(1:3);

            % Since Tmat is a known function of variable, use the analytical inverse
            % of matPFinal to avoid numerical issues:
            invmatPFinal = [10/obj.Variable^3 -4/obj.Variable^2 1/(2*obj.Variable);
                            -15/obj.Variable^4 7/obj.Variable^3 -1/obj.Variable^2;
                            6/obj.Variable^5 -3/obj.Variable^4, 1/(2*obj.Variable^3)];

            coeffVec(4:6) = invmatPFinal * B;

            coefficients = coeffVec';
        end

        function coefficients = computeQuarticTrajectory(obj)
        %computeTrajectory Computes the coefficients of the polynomial
        %   This function computes the quartic polynomial coefficients
        %   given the zeroth, first, and second derivative at P = 0 and
        %   P = finalVariable

        % First three coefficients can be found by setting P=0 using
        % boundary conditions
            coeffVec = [obj.StartCondition(1) obj.StartCondition(2) obj.StartCondition(3)/2.0 0 0]';

            % The last three coefficients may be found by substituting the first three
            % into the expressions for zeroth, first, and second derivative at
            % P=finalVariable:
            % [fd(P) fdd(P)]' = matP0*[C1 C2 C3]' + matPFinal*[C4 C5]
            % where matPFinal =
            %   [3*finalVariable^2  4*finalVariable^3; ...
            %    6*finalVariable    12*finalVariable^2];
            % This expression may be solve for C4, and C5

            matP0 = [0 1 2*obj.Variable; 0 0 2];

            B = [obj.EndCondition(2); obj.EndCondition(3)] - matP0*coeffVec(1:3);

            % Since Tmat is a known function of variable, use the analytical inverse
            % of matPFinal to avoid numerical issues:
            invmatPFinal = 1/(12*obj.Variable^4) * [12*obj.Variable^2 -4*obj.Variable^3; ...
                                -6*obj.Variable 3*obj.Variable^2];

            coeffVec(4:5) = invmatPFinal * B;

            coefficients = [coeffVec', 0];
        end
    end
end
