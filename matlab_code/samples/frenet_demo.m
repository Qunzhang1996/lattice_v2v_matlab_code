% frenet_demo.m

clc;
waypoints = [0 0 -0.3531; ...
            50 -10 0.1148; ...
            100 10 0.4193; ...
            150 20 -0.0685
            200 10 -0.1484
            250 10 0.0742];
refPath = FrenetReferencePath(waypoints(:,1:3));

ss = 0:250;
p1s = zeros(length(ss),2);
p2s = zeros(length(ss),2);
for i = 1:length(ss)
    p1 = refPath.frenet2global([ss(i), 0, 0,2, 0, 0]);
    p2 = refPath.frenet2global([ss(i), 0, 0,-2, 0, 0]);
    p1s(i, :) = p1(1:2);
    p2s(i, :) = p2(1:2);
end
xy = [52, 0];
fstate = refPath.global2frenet([xy, 0, 0, 0, 0]);
rp = refPath.interpolate(fstate(1));
pp = refPath.interpolate(ss);
hold on;
plot(xy(1),xy(2), 'ro');
plot(rp(1),rp(2), 'bo');
plot([xy(1),rp(1)],[xy(2), rp(2)], '--');
plot(pp(:,1),pp(:,2), 'color','r')
plot(waypoints(:,1),waypoints(:,2), 'o','MarkerFaceColor','r')
plot(p1s(:,1),p1s(:,2), 'marker','.', 'color','k')
plot(p2s(:,1),p2s(:,2), 'marker','.', 'color','k')
axis equal;