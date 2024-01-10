classdef TrajPredicter < handle

    properties(Access=public)
        numActors               % 演员车数量
        filters                 % 滤波器
        dT                      % 时间间隔
        refPath                 % 道路参考线
    end
    methods
        function obj = TrajPredicter(curActorState, dT, refPath)
            obj.numActors = size(curActorState,1);
            obj.filters   = cell(1,obj.numActors);
            obj.dT        = dT;
            obj.refPath   = refPath;
            for i = 1:obj.numActors
                x     = curActorState(i,1);
                y     = curActorState(i,2);
                theta = curActorState(i,3);
                speed = curActorState(i,4);
                measurement = [x, speed*cos(theta), y, speed*sin(theta)];
                obj.filters{i}   = InitFrenetStateFilters(measurement.', refPath);
            end
        end
        function update(obj, curActorState, refPath)
            for i = 1:obj.numActors
                x     = curActorState(i,1);
                y     = curActorState(i,2);
                theta = curActorState(i,3);
                speed = curActorState(i,4);
                measurement = [x, speed*cos(theta), y, speed*sin(theta)];
                [xpred, ~] = obj.filters{i}.predict(obj.dT);%#ok
                [xest, ~]  = obj.filters{i}.correct(measurement.', refPath);%#ok
            end
        end
        function futureTrajectory = getFutureTrajectory(obj)
            futureTrajectory = repelem(struct('Trajectory',[]),obj.numActors,1);
            for i = 1:obj.numActors
                futureTrajectory(i).Trajectory = zeros(30,6);
                filter = obj.filters{i}.clone();
                for j = 1:30
                    [xpred,~] = filter.predict(obj.dT);
                    globalstate = filterToGlobalState(obj, xpred);
                    futureTrajectory(i).Trajectory(j,:) = globalstate.';
                end
            end
        end

        function globalState = filterToGlobalState(obj,filterState)

            % Assemble as Frenet state to use frenet2global
            numStates = size(filterState,2);
            s = filterState(1,:);
            ds = filterState(2,:);
            dds = zeros(1,numStates);
            d = filterState(3,:);
            dd = filterState(4,:);
            ddbyds = zeros(1,numStates);
            ddbyds(abs(ds) > 0.5) = dd(abs(ds) > 0.5)./ds(abs(ds) > 0.5);
            dd2ds2 = zeros(1,numStates);
            frenetState = [s;ds;dds;d;ddbyds;dd2ds2];

            % Convert to global state
            % refPath = helperGetReferencePath;
            globalState = frenet2global(obj.refPath,frenetState')';
        end
    end
end
