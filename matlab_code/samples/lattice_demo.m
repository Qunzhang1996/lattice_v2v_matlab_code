%lattice_demo.m

clc;
waypoints = [0 0 -0.3531; ...
            50 -10 0.1148; ...
            100 10 0.4193; ...
            150 20 -0.0685
            200 10 -0.1484
            250 10 0.0742];
refPath = FrenetReferencePath(waypoints(:,1:3));
pp = refPath.interpolate(0:0.1:refPath.Length);
hold on;
plot(pp(:,1),pp(:,2), 'color','r')
connector  = TrajectoryGeneratorFrenet(refPath);
for i = -3:3
    initState = [0 6 0 0 0 0];  % [S ds ddS L dL ddL]
    termState = [NaN 6 0 i 0 0]; % [S ds ddS L dL ddL]
    [~,trajGlobal] = connect(connector,initState,termState,5);
    plot(trajGlobal.Trajectory(:,1),trajGlobal.Trajectory(:,2), 'marker', '.');
end
axis equal;
