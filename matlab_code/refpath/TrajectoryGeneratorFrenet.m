classdef TrajectoryGeneratorFrenet < handle
    properties
        %ReferencePath Reference path for generated trajectories
        %
        %   A referencePathFrenet object used to map trajectories generated
        %   between Frenet states into the global coordinate system.
        ReferencePath
        DefaultTimeResolution = 0.1;
        %TimeResolution Discretization time interval
        TimeResolution = 0.1;
        end
    methods
        function obj = TrajectoryGeneratorFrenet(referencePath, varargin)
        %trajectoryGeneratorFrenet Generates alternate trajectories relative to a reference path
            obj.ReferencePath = referencePath;
            
            % Parse optional arguments
            if nargin > 1
                validatestring(varargin{1},{'TimeResolution'},'trajectoryGeneratorFrenet');
                
                % Set optional properties
                obj.TimeResolution = varargin{2};
            end
        end
        
        function [frenetTrajectory, globalTrajectory] = connect(obj, initialState, terminalState, timeSpan)
        %connect Connect initial and terminal Frenet states
        %
        %   FRENETTRAJECTORY = connect(obj, INITIALSTATE, TERMINALSTATE, TIMESPAN)
        %   connects initial states to terminal states over a span of time.
        %   INITIALSTATE and TERMINALSTATE are 6-column matrices of Frenet 
        %   states, [S dS ddS L dL ddL], and TIMESPAN is a scalar or vector
        %   of positive time spans over which the initial and terminal states
        %   are connected. 
        %
        %   This function supports generation of 1-to-N, N-to-1, and N-to-N
        %   (pairwise) trajectories based on the number of states and times.
        %   FRENETTRAJECTORY is an Nx1 struct array. Each element contains
        %   Trajectory, an Mx6 matrix of Frenet states, and Time, an Mx1
        %   vector of times at which each states occurs.
        %
        %   To calculate the longitudinal terminal state based on a 
        %   4th-order polynomial while satisfying all other boundary 
        %   conditions, specify the arc length, S, as NaN.
        %   
        %   [__,GLOBALTRAJECTORY] = connect(obj, INITIALSTATE, TERMINALSTATE, TIMESPAN)
        %   returns an optional output, GLOBALTRAJECTORY, an Nx1 struct array. 
        %   Each element contains Time, defined above, and Trajectory, an
        %   Mx6 matrix of global states, [x,y,theta,kappa,speed,acceleration].
        %
        %   Example:
        %         % Generate a reference path from a set of waypoints.
        %         refPath = referencePathFrenet([0 0; 50 20; 0 50; -10 0]);
        % 
        %         % Create a trajectory generator object.
        %         connector = trajectoryGeneratorFrenet(refPath);
        % 
        %         % Generate a 5 second trajectory between the path origin
        %         % and a point 30m down the path.
        %         initState = [ 0 0 0 0 0 0]; % [S ds ddS L dL ddL]
        %         termState = [30 0 0 0 0 0]; % [S ds ddS L dL ddL]
        %         trajFrenet = connect(connector, initState, termState, 5);
        % 
        %         % Generate trajectories that traverse the same arc length, 
        %         % but deviate laterally from the reference path.
        %         termStateDeviated = termState + (-3:3)'*[0 0 0 1 0 0];
        %         [trajFrenet, trajGlobal] = connect(connector, initState, termStateDeviated, 5);
        % 
        %         % Specify a terminal state with 5m arc length, 10m/s 
        %         % velocity, and 10m longitudinal offset. These boundary 
        %         % conditions are not realistic, and need adjustment in 
        %         % the next steps.
        %         initState = [0 0 0 0 0 0];
        %         unrealTermState = [5 10 0 10 0 0];
        %         [~, unrealTrajGlobal] = connect(connector, initState, unrealTermState, 3);
        %
        %         % Display the trajectory and note how the vehicle backs
        %         % up before accelerating to the terminal state. 
        %         show(refPath);
        %         hold on;
        %         axis equal;
        %         plot(unrealTrajGlobal.Trajectory(:,1),unrealTrajGlobal.Trajectory(:,2),'b');
        % 
        %         % Remove the terminal longitudinal constraint.
        %         relaxedTermState = [NaN 10 0 10 0 0];
        %
        %         % Generate a new trajectory.
        %         [~, trajGlobalRelaxed] = connect(connector, initState, relaxedTermState, 3);
        %         plot(trajGlobalRelaxed.Trajectory(:,1),trajGlobalRelaxed.Trajectory(:,2),'g');
        %         legend({'Waypoints','ReferencePath','Unrealistic Trajectory', ...
        %           'Relaxed 4th order trajectory'})
        %
        %       See also referencePathFrenet, trajectoryGeneratorFrenet
            
            narginchk(4,4);
            
            % Validate inputs
            obj.validateConnect(initialState, terminalState, timeSpan);
            
            % Organize inputs into two sets of N-row state pairs
            [f0,f1,tF] = obj.parseTrajInputs(initialState, terminalState, timeSpan);
            
            % Generate trajectories in Frenet coordinates
            frenetTrajectory = obj.connectPairs(f0, f1, obj.TimeResolution, tF);
            
            if nargout == 2
                % Convert trajectories to global coordinates if requested
                globalTrajectory = frenetTrajectory;
                for i = 1:numel(frenetTrajectory)
                    globalTrajectory(i).Trajectory = obj.ReferencePath.frenet2global(frenetTrajectory(i).Trajectory);
                end
            end
        end
        function set.ReferencePath(obj, refPathObj)
        %set.ReferencePath
            obj.ReferencePath = refPathObj;
        end
        
        function set.TimeResolution(obj, resolution)
        %set.TimeResolution
            obj.TimeResolution = resolution;
        end
    end
    
    methods (Access = public, Static)
        function trajectories = connectPairs(f0,f1,dt,timespan)
        %connectPairs Connect N pairs of start/end states
            % Number of expected trajectories
            numTraj = size(f0,1);
            
            % Preallocate trajectory variables
            longitudinalTrajectories = repmat({QuinticQuarticTrajectory([0 0 0],[1 1 1],1)},numTraj,1);
            lateralTrajectories = repmat({QuinticQuarticTrajectory([0 0 0],[1 1 1],1)},numTraj,1);
            
            traj = [];
            t = [];
            coder.varsize('traj',[inf, 6],[1 0]);
            coder.varsize('t',[inf,1],[1 0]);
            
            trajectories = repmat(struct('Trajectory',traj,'Times',t), numTraj, 1);
            coder.varsize('trajectories',[inf 1],[1 0]);
            for trajIdx = 1:numTraj
                % Extract longitudinal boundary conditions
                sV0 = f0(trajIdx,1:3);
                sV1 = f1(trajIdx,1:3);
                tFinal = timespan(trajIdx);

                % Fit quintic/quartic polynomial
                longitudinalTrajectories{trajIdx} = ...
                    QuinticQuarticTrajectory(sV0,sV1,tFinal);

                % Evaluate longitudinal trajectory at evenly spaced
                % intervals wrt time
                t = (dt:dt:dt*(ceil(tFinal/dt)+1))'-dt;
                sV = longitudinalTrajectories{trajIdx}.evaluate([0 1 2], t)';

                % Extract longitudinal boundary conditions
                dV0 = f0(trajIdx,4:6);
                dV1 = f1(trajIdx,4:6);
                dsMax = sV(end,1)-sV(1);
                
                % Fit quintic/quartic polynomial
                lateralTrajectories{trajIdx} = ...
                    QuinticQuarticTrajectory(dV0,dV1,dsMax);

                % Evaluate lateral trajectory based on arclength
                ds = sV(:,1)-sV(1);
                dV = lateralTrajectories{trajIdx}.evaluate([0 1 2], ds)';

                % Combine and populate output struct
                trajectories(trajIdx).Trajectory = reshape([sV dV],[],6);
                trajectories(trajIdx).Times = t;
            end
        end
        
        function [x0,x1,tF] = parseTrajInputs(initState, termState, time)
        %parseTrajInputs Parses and organizes the inputs into two sets of N-row state pairs
        
            % Verify the number of state/time inputs is supported
            numInit = size(initState,1);
            numTerm = size(termState,1);
            numTime = numel(time);
            x0 = initState;
            x1 = termState;
            tF = time(:);
            
            numTraj = max([numInit, numTerm, numTime]);
            if numTraj == 1 || all(numTraj == [numInit numTerm numTime])
                return;
            else
                if numInit ~= numTraj
                    % Expand initial state
                    coder.internal.errorIf(numInit ~= 1, 'nav:navalgs:trajectorygeneratorfrenet:MismatchedNumStates');
                    coder.varsize('x0',[inf 6],[1 0]);
                    x0 = repmat(x0, numTraj, 1);
                end
                if numTerm ~= numTraj
                    % Expand terminal state
                    coder.internal.errorIf(numTerm ~= 1, 'nav:navalgs:trajectorygeneratorfrenet:MismatchedNumStates');
                    coder.varsize('x1',[inf 6],[1 0]);
                    x1 = repmat(x1, numTraj, 1);
                end
                if numTime ~= numTraj
                    % Expand time
                    coder.internal.errorIf(numTime ~= 1, 'nav:navalgs:trajectorygeneratorfrenet:MismatchedNumTimes');
                    coder.varsize('tF',[inf 1],[1 0]);
                    tF = repmat(tF, numTraj, 1);
                end
            end
        end
        
        function validateConnect(initState, termState, time)
        %validateConnect Validate type/size of inputs to connect method
        %
        %   The first column of terminal state is allowed to be nan, 
        %   reducing boundary constraints by 1 and allowing us to use 4th
        %   order polynomial to satisfy velocity, acceleration, and initial
        %   position.
            
            % Verify attributes of initial state
            validateattributes(initState,{'numeric'},{'size',[nan 6],'finite','nonnan'},'connect','initState');
            
            % Verify size of terminal states
            assert(size(termState,2) == 6);
            % Verify terminal states form valid 4th/5th order boundary conditions
            validateattributes(termState(:,2:end),{'numeric'},{'nonnan','finite'},'connect','termState');
            
            % Verify the times are valid
            validateattributes(time,{'numeric'},{'finite','nonempty','vector','positive'},'connect','time');
        end
    end
end
