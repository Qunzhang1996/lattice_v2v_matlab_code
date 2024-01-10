classdef ScenarioEnv < handle
    
    properties(Access=public)
        actors                  % actors信息
        refPath                 % 参考线
        laneNum                 % 车道数量
        laneCenters             % 车道中心点偏移
        laneBounds              % 车道边界
        laneWidth               % 车道宽
        egoInitState            % ego车初始状态
        numActors               % 演员车数量
        egop                    % ego绘图句柄
        carp                    % actors绘图句柄
        cartrajp                % actors轨迹绘图句柄
        trajp                   % 轨迹绘图句柄
        trajopt                 % 最优轨迹绘图句柄
        SampleTime              % 采样时间
        replanRate              % 重规划频率
        maxHorizon              % 滚动窗口长度
        fig                     % 图窗figure句柄
        ax                      % 图窗axes句柄
        curState                % actors的当前状态
    end
    methods
        function updateShow(obj, egoState, curActorState, globalTraj,optimalTrajectory, futureTrajectory)% 更新绘图
            curpos = egoState(1:3);
            xy = [-2.5 -2.5 2.5 2.5
                   1   -1  -1   1];
            M = [cos(curpos(3)), -sin(curpos(3));sin(curpos(3)) cos(curpos(3))];
            xy = M*xy+curpos(1:2)';
            set(obj.egop,'xdata',xy(1,:),'ydata',xy(2,:));
            for idx = 1:obj.numActors
                curpos = curActorState(idx,1:3);
                xy = [-2.5 -2.5 2.5 2.5
                    1   -1  -1   1];
                M = [cos(curpos(3)), -sin(curpos(3));sin(curpos(3)) cos(curpos(3))];
                xy = M*xy+curpos(1:2)';
                set(obj.carp(idx),'xdata',xy(1,:),'ydata',xy(2,:));
            end

            for idx = 1:12
                if idx > length(globalTraj)
                    set(obj.trajp(idx), 'xdata',nan, 'ydata',nan);
                else
                    set(obj.trajp(idx), 'xdata',globalTraj(idx).Trajectory(:,1), 'ydata',globalTraj(idx).Trajectory(:,2));
                end
            end
            set(obj.trajopt, 'xdata',optimalTrajectory(:,1), 'ydata',optimalTrajectory(:,2));

            for idx = 1:obj.numActors
                xs = futureTrajectory(idx).Trajectory(:,1);
                ys = futureTrajectory(idx).Trajectory(:,2);
                set(obj.cartrajp(idx),'xdata',xs,'ydata',ys);
            end
        end
        function curState = getActorInfo(obj)% 获取actors的信息
            curState = obj.curState;
        end
        function update(obj) %actors位置和未来轨迹更新
            % 获取非ego的位姿和未来轨迹
            numActor  = obj.numActors;
            curState1 = zeros(numActor,6);
            for k = 1:numActor
                obj.actors{k}.s = obj.actors{k}.s + obj.SampleTime*obj.actors{k}.speed(1);
            end
            poses1 = struct('Position',[],'Velocity',[],'Yaw',[],'AngularVelocity',[]);
            poses  = repmat(poses1, 1, numActor);
            for k = 1:numActor
                p1 = obj.actors{k}.refPath.interpolate(obj.actors{k}.s);
                poses1.Position = [p1(1:2),0];
                poses1.Velocity = obj.actors{k}.speed(1)*[cos(p1(3)),sin(p1(3)),0];
                poses1.Yaw = p1(3)*180/pi;
                poses1.AngularVelocity = [0,0,p1(4)*obj.actors{k}.speed(1)*180/pi];
                poses(k) = poses1;
            end
            for j = 1:numActor
                actIdx = j;
                xy = poses(actIdx).Position(1:2);
                v  = norm(poses(actIdx).Velocity,2);
                th = atan2(poses(actIdx).Velocity(2),poses(actIdx).Velocity(1));
                k = poses(actIdx).AngularVelocity(3)/v/180*pi;
                curState1(j,:) = [xy th k v 0];
            end
            obj.curState = curState1;
        end
        function show(obj) % 绘图
            % 绘制场景
            hold on;
            ss = 0:1:obj.refPath.Length;
            len = length(ss);
            for i = 1:length(obj.laneCenters)
                laneCenter = obj.laneCenters(i);
                frestate = [0, 0, 0, laneCenter-obj.laneWidth/2, 0, 0];
                frestates = repmat(frestate, len, 1);
                frestates(:,1) = ss';
                pp1 = obj.refPath.frenet2global(frestates);
                plot(pp1(:,1),pp1(:,2), 'k','LineWidth',1,'LineStyle','--');

                frestate = [0, 0, 0, laneCenter+obj.laneWidth/2, 0, 0];
                frestates = repmat(frestate, len, 1);
                frestates(:,1) = ss';
                pp1 = obj.refPath.frenet2global(frestates);
                plot(pp1(:,1),pp1(:,2), 'k','LineWidth',1,'LineStyle','--');
            end
            laneCenter = obj.laneBounds(1);
            frestate = [0, 0, 0, laneCenter, 0, 0];
            frestates = repmat(frestate, len, 1);
            frestates(:,1) = ss';
            pp1 = obj.refPath.frenet2global(frestates);
            plot(pp1(:,1),pp1(:,2), 'k','LineWidth',2,'LineStyle','-');
            laneCenter = obj.laneBounds(end);
            frestate = [0, 0, 0, laneCenter, 0, 0];
            frestates = repmat(frestate, len, 1);
            frestates(:,1) = ss';
            pp1 = obj.refPath.frenet2global(frestates);
            plot(pp1(:,1),pp1(:,2), 'k','LineWidth',2,'LineStyle','-');
            axis equal;

            obj.fig = gcf;
     
            obj.ax = gca;
            obj.egop = patch(nan,nan,'r');
            obj.carp = zeros(obj.numActors,1);
            obj.cartrajp = zeros(obj.numActors,1);
            %绘制演员车真实轨迹
            for i = 1:obj.numActors
                obj.carp(i) = patch(nan,nan,'g');
                obj.cartrajp(i) = plot(nan,nan,'k','LineWidth',1,'LineStyle','-.');
                % endp = obj.actors{i}.refPath.interpolate(realmax);
                % wss  = obj.actors{i}.refPath.interpolate(0:1:endp(6));
                % plot(wss(:,1),wss(:,2),'k');
            end
            obj.trajp = zeros(12,1);
            for i = 1:12
                obj.trajp(i) = plot(nan,nan,'b');
            end
            title('Lattice Planner Demo');
            % 绘制演员车预测轨迹
            obj.trajopt = plot(nan,nan, 'marker','.', 'color','g','linewidth',1.5);
            set(gcf, 'WindowState','maximized');
        end
        function [laneCenters,laneBounds]= calcLaneCenters(obj)
            % 计算车道中心点的偏移
            laneCenters = zeros(1,obj.laneNum);
            laneBounds  = zeros(1,obj.laneNum+1);
            if mod(obj.laneNum,2) == 0
                % 偶数车道
                for i = 1:obj.laneNum
                    laneCenters(i) = (obj.laneNum/2 - 0.5)*obj.laneWidth - (i-1)*obj.laneWidth;
                    laneBounds(i)  = laneCenters(i) + 0.5*obj.laneWidth;
                end
            else
                % 奇数车道
                for i = 1:obj.laneNum
                    laneCenters(i) = (obj.laneNum-1)/2*obj.laneWidth - (i-1)*obj.laneWidth;
                    laneBounds(i)  = laneCenters(i) + 0.5*obj.laneWidth;
                end
            end
            laneBounds(end) = laneBounds(end-1) - obj.laneWidth;
        end
        function loadScenarioFromFile(obj, filename)
            data = load(filename);
            % 加载道路中心坐标
            waypoints = data.data.RoadSpecifications.Centers(:,1:2);
            % 加载车道数量
            obj.laneNum     = data.data.RoadSpecifications.Lanes.NumLanes;
            % 车道宽度，目前仅支持车道宽度一致
            Widths = data.data.RoadSpecifications.Lanes.Width;
            obj.laneWidth   = Widths(1);
            assert(all(abs(Widths-Widths(1)) < 0.1),'目前仅支持车道宽度一致');
            [obj.laneCenters,obj.laneBounds] = calcLaneCenters(obj);
            obj.refPath     = FrenetReferencePath(waypoints);
            % x y theta kappa speed a
            Actors = data.data.ActorSpecifications;
            egoid  = data.data.EgoCarId;
            ActorIDs = arrayfun(@(x) x.ActorID, Actors);
            egoidx   = find(ActorIDs == egoid);
            assert(length(egoidx) == 1,'egoid只能有1个');
            obj.egoInitState    = [Actors(egoidx).Position,0,0,0];
            % 添加演员车，每个演员车需要有初始位置，路点，速度
            cnt = 0;
            obj.actors = cell(1,length(Actors)-1);
            for i = 1:length(Actors)
                if Actors(i).ActorID ~= egoid
                    car1.Position = Actors(i).Position;
                    car1.waypoints = Actors(i).Waypoints;
                    car1.speed = Actors(i).Speed;
                    car1.refPath = FrenetReferencePath(car1.waypoints(:,1:2));
                    car1.s = 0;
                    cnt = cnt + 1;
                    obj.actors{cnt} = car1;
                end
            end
            obj.numActors = length(obj.actors);
            for j = 1:obj.numActors
                xyth = obj.actors{j}.Position;
                v  = obj.actors{j}.speed(1);
                k = 0;
                obj.curState(j,:) = [xyth k v 0];
            end
        end
        function loadScenarioFromCode(obj)

            % 需要定义道路中心
            waypoints = [0 50; 150 50; 300 75; 310 75; 400 0; 300 -50; 290 -50; 0 -50]; % in meters
            % 车道数量
            obj.laneNum     = 4;
            % 车道宽度，目前仅支持车道宽度一致
            obj.laneWidth   = 3.6;
            [obj.laneCenters,obj.laneBounds] = calcLaneCenters(obj);
            obj.refPath     = FrenetReferencePath(waypoints);
            % x y theta kappa speed a
            obj.egoInitState    = [0, 48.204, 0, 0, 0, 0];
            obj.egoInitState    =[0.1211   51.7959   -0.0673         0         0         0];
            % 添加演员车，每个演员车需要有初始位置，路点，速度
            car1.Position = [34.7 49.3 0];
            car1.waypoints = [34.7 49.3 0;
                60.1 48.2 0;
                84.2 47.9 0;
                119 49.3 0;
                148.1 51.4 0;
                189.6 58.7 0;
                230.6 68 0;
                272.6 74.7 0;
                301.4 77.5 0;
                316.7 76.8 0;
                332.4 75.2 0;
                348.9 72.2 0;
                366.2 65.1 0;
                379.6 55.6 0];

            car1.speed = [10;10;10;10;10;10;10;10;10;10;10;10;10;10];

            car2.Position = [17.6 46.7 0];
            car2.waypoints = [17.6 46.7 0;
                43.4 45.5 0;
                71.3 43.8 0;
                102.3 43.5 0;
                123.5 45.5 0;
                143.6 47.4 0;
                162.4 50 0;
                198.5 61 0;
                241.1 70.1 0;
                272.3 74.1 0;
                292 76.6 0;
                312.8 77.2 0;
                350.3 75.2 0;
                362.5 70.4 0;
                375.9 63.3 0;
                390.7 49.9 0;
                401.3 33 0];
            car2.speed = [9;9;9;9;9;9;9;9;9;9;9;9;9;9;9;9;9];

            car3.Position = [62.6 51.9 0];
            car3.waypoints = [62.6 51.9 0;
                87.4 51.3 0;
                117.7 52.2 0;
                147.6 55 0;
                174.9 59.7 0;
                203.3 65.8 0;
                265 69.7 0;
                288.3 73.1 0;
                314.5 73.1 0;
                334.9 70.8 0;
                360 59.9 0];
            car3.speed = [6;6;6;6;6;6;6;6;6;6;6];

            car4.Position = [101.7 41.1 0];
            car4.waypoints = [101.7 41.1 0;
                124.6 42 0;
                148.5 43.9 0;
                171.9 48.2 0;
                197.1 52.8 0;
                222.3 58.5 0;
                252.4 64.4 0;
                281.4 68.5 0;
                307.7 69.5 0;
                329.9 68.2 0;
                352.7 62.8 0];
            car4.speed = [7;7;7;7;7;7;7;7;7;7;7];

            car5.Position = [251.3 75.6 0];
            car5.waypoints = [251.3 75.6 0;
                255.7 76.7 0];
            car5.speed = [0.01;0.01];
            car1.refPath = FrenetReferencePath(car1.waypoints(:,1:2));
            car2.refPath = FrenetReferencePath(car2.waypoints(:,1:2));
            car3.refPath = FrenetReferencePath(car3.waypoints(:,1:2));
            car4.refPath = FrenetReferencePath(car4.waypoints(:,1:2));
            car5.refPath = FrenetReferencePath(car5.waypoints(:,1:2));
            car1.s = 0;
            car2.s = 0;
            car3.s = 0;
            car4.s = 0;
            car5.s = 0;
            obj.actors = {car1,car2,car3,car4,car5};
            obj.numActors = length(obj.actors);
            for j = 1:obj.numActors
                xyth = obj.actors{j}.Position;
                v  = obj.actors{j}.speed(1);
                k = 0;
                obj.curState(j,:) = [xyth k v 0];
            end
            
        end
        function obj = ScenarioEnv()% 构造函数
        end
    end


end