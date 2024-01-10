clear;
figure
waypoints = [0 50; 150 50; 300 75; 310 75; 400 0; 300 -50; 290 -50; 0 -50]; % in meters

% car1.Position = [34.7 49.3 0];
% car1.s = 0;
% car1.waypoints = [34.7 49.3 0;
%     60.1 48.2 0;
%     84.2 47.9 0;
%     119 49.3 0;
%     148.1 51.4 0;
%     189.6 58.7 0;
%     230.6 68 0;
%     272.6 74.7 0;
%     301.4 77.5 0;
%     316.7 76.8 0;
%     332.4 75.2 0;
%     348.9 72.2 0;
%     366.2 65.1 0;
%     379.6 55.6 0];
% car1.speed = [10;10;10;10;10;10;10;10;10;10;10;10;10;10];
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
car1 = car2;
car1.s =  0;
hold on;
car1p = patch(nan,nan,'g');
estp  = plot(nan,nan,'ro');
refPath  = referencePathFrenet(waypoints);
ss       = refPath.interpolate(0:1:refPath.PathLength);
refPath1 = referencePathFrenet(car1.waypoints(:,1:2));
ss1      = refPath1.interpolate(0:1:refPath1.PathLength);
plen = refPath.PathLength;
pss  = (0:plen)';
gss  = refPath.frenet2global([pss, zeros(length(pss),2),zeros(length(pss),1)+1.8,zeros(length(pss),2)]);
plot(gss(:,1),gss(:,2),'k--');
gss  = refPath.frenet2global([pss, zeros(length(pss),2),zeros(length(pss),1)-1.8,zeros(length(pss),2)]);
plot(gss(:,1),gss(:,2),'k--');
plot(ss(:,1),ss(:,2),'k');
plot(ss1(:,1),ss1(:,2));
axis equal;
dT = 0.1;
curActorState = refPath1.interpolate(car1.s);
x     = curActorState(1);
y     = curActorState(2);
theta = curActorState(3);
speed = car1.speed(1);
measurement1 = [x, speed*cos(theta), y, speed*sin(theta)];
filter   = InitFrenetStateFilters(measurement1.', refPath);
State1 = filter.State;
pxys = zeros(30,2)*nan;
trajp = plot(pxys(:,1),pxys(:,2),'k','Marker','.');
for i = 1:5
    preh(i) = plot(nan,nan,'r','Marker','.');
end
% plot(car1.waypoints(:,1),car1.waypoints(:,2),'k')
data = load("pre_data.mat");
for i = 1:168
    % 更新car
    car1.s = car1.s + car1.speed(1) * dT;
    p1 = refPath1.interpolate(car1.s);
    xy = p1(1:2);
    v  = car1.speed(1);
    th = p1(3);
    k = p1(4);
    curState1 = [xy th k v 0];
    curpos = curState1(1:3);
    xy = [-2.5 -2.5 2.5 2.5
        1   -1  -1   1];
    M = [cos(curpos(3)), -sin(curpos(3));sin(curpos(3)) cos(curpos(3))];
    xy = M*xy+curpos(1:2)';
    set(car1p,'xdata',xy(1,:),'ydata',xy(2,:));

    x     = curState1(1);
    y     = curState1(2);
    theta = curState1(3);
    speed = curState1(5);
    measurement = [x, speed*cos(theta), y, speed*sin(theta)];
    [xpred, ~] = filter.predict(dT);
    [xest, ~]  = filter.correct(measurement.', refPath);
    globalState=filterState2global(xest,refPath);
    set(estp, 'xdata',globalState(1), 'ydata',globalState(2));
    tmpfilter = filter.clone();
    for j = 1:30
        [xpred, ~] = tmpfilter.predict(dT);
        globalState=filterState2global(xpred,refPath);
        pxys(j,:) = globalState(1:2);
    end
    set(trajp, 'xdata',pxys(:,1), 'ydata',pxys(:,2));

    tgts = data.pre_data{i}.tgtPoses;
    for idx = 1:length(tgts)
        set(preh(idx), 'xdata',tgts(idx).States(:,1),'ydata',tgts(idx).States(:,2));
    end
    drawnow;
end

function  globalState=filterState2global(xest,refPath)
    s = xest(1);
    ds = xest(2);
    dds = 0;
    d = xest(3);
    dd = xest(4);
    ddbyds = 0;
    ddbyds(abs(ds) > 0.5) = dd(abs(ds) > 0.5)./ds(abs(ds) > 0.5);
    dd2ds2 = 0;
    frenetState = [s;ds;dds;d;ddbyds;dd2ds2];
    globalState = frenet2global(refPath,frenetState')';
end



