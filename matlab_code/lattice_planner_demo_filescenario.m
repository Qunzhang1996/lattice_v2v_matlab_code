% lattice_planner_demo.m

clear;
% 车辆参数
carLen   = 4.7; 
carWidth = 1.8; 
carRear = 1.175; 
% 场景
scenario = ScenarioEnv;
% 加载代码中定义的场景并绘图
scenario.loadScenarioFromFile('scenario1.mat');
scenario.show;
% 获取车道宽度，车道边界，车道中心
laneWidth   = scenario.laneWidth;
laneBounds  = scenario.laneBounds;
laneCenters = scenario.laneCenters;
% 获取参考线
refPath   = scenario.refPath;
connector = TrajectoryGeneratorFrenet(refPath);
scenario.SampleTime = connector.TimeResolution;

% 动态胶囊参数
Geometry.Length = carLen; % in meters
Geometry.Radius = carWidth/2; % in meters
Geometry.FixedTransform = -carRear; % in meters
% 规划参数
replanRate = 10; % Hz
scenario.replanRate = replanRate;
timeHorizons = 1:3; % in seconds
maxHorizon = max(timeHorizons); % in seconds
scenario.maxHorizon = maxHorizon;
latDevWeight    =  1;
timeWeight      = -1;
speedWeight     =  1;
maxAcceleration =  15; % in meters/second^2
maxCurvature    =   1; % 1/meters, or radians/meter
minVelocity     =   0; % in meters/second
speedLimit = 11; % in meters/second
safetyGap  = 10; % in meters

% 初始化
numActors = scenario.numActors;
% 非ego车辆当前位姿
actorPoses = repelem(struct('States',[]),numActors,1);
% ego起点
egoState = scenario.egoInitState;
curActorState = scenario.getActorInfo();
trajPredicter = TrajPredicter(curActorState, scenario.SampleTime, refPath);
% 打开视频准备保存
outmp4 = VideoWriter('lattice_demo.mp4','MPEG-4');
open(outmp4);
% 主循环
for loop = 1:1500
    % 场景更新
    scenario.update();
    % 获取场景内Actor的信息
    curActorState = scenario.getActorInfo();
    % use actor's [x,y,th,v] to predict future trajectory
    trajPredicter.update(curActorState(:,[1,2,3,5]), refPath);
    futureTrajectory = trajPredicter.getFutureTrajectory();
    % 轨迹终点状态采样
    [allTS, allDT, numTS] = SamplingEndcontions(refPath, laneWidth, laneBounds, egoState, curActorState, safetyGap, speedLimit, timeHorizons);
    % 轨迹评估
    costTS = EvaluateTSCost(allTS,allDT,laneWidth,laneCenters,speedLimit,speedWeight, latDevWeight, timeWeight);
    % global转frenet
    egoFrenetState = refPath.global2frenet(egoState);
    % 轨迹生成
    [frenetTraj,globalTraj] = connector.connect(egoFrenetState,allTS,allDT);
    % 轨迹筛选
    isValid = EvaluateTrajectory(globalTraj,maxAcceleration,maxCurvature,minVelocity);
    % 碰撞检测--设置actor的未来轨迹
    for i = 1:numActors
        actorPoses(i).States = futureTrajectory(i).Trajectory(:,1:3);
    end
    % 按照代价对轨迹进行排序
    [cost, idx] = sort(costTS);
    % 开始轨迹
    optimalTrajectory = [];
    for i = 1:numel(idx)
        % 根据筛选的结果选取有效的轨迹
        if isValid(idx(i))
            % 设置ego的未来轨迹
            egoPoses.States = globalTraj(idx(i)).Trajectory(:,1:3);
            % 碰撞检测
            isColliding = checkTrajCollision(egoPoses, actorPoses, Geometry);
            if all(~isColliding)
                % 无碰撞，则找到最优的全局路径
                optimalTrajectory = globalTraj(idx(i)).Trajectory;
                break;
            end
        end
    end
    if isempty(optimalTrajectory)
        % 如果没有找到轨迹，则报错
        error('No valid trajectory has been found.');
    else
        % 更新绘图
        scenario.updateShow(egoState, curActorState,globalTraj, optimalTrajectory, futureTrajectory);
        % 更新ego位置
        egoState = optimalTrajectory(2,:);
    end
    % 立即绘图
    drawnow;
    % 延时
    pause(0.1);
    frame = getframe(gcf);
    writeVideo(outmp4,frame);
    writeVideo(outmp4,frame);
    writeVideo(outmp4,frame);
    writeVideo(outmp4,frame);
end
close(outmp4);