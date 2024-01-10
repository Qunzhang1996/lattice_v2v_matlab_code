% refPath_demo.m

clc;
waypoints = [0   0 -0.3531; ...
            50 -10 0.1148; ...
            100 10 0.4193; ...
            150 20 -0.0685
            200 10 -0.1484
            250 10 0.0742];
refPath = FrenetReferencePath(waypoints(:,1:3));
ss = 0:260;
pp = refPath.interpolate(ss);
hold on;
plot(waypoints(:,1),waypoints(:,2), 'ro');
plot(pp(:,1),pp(:,2), 'color','r')

figure
subplot(3,1,1)
plot(pp(:,6), pp(:,3))
title('航向角')
subplot(3,1,2)
plot(pp(:,6), pp(:,4))
title('曲率')
subplot(3,1,3)
plot(pp(:,6), pp(:,5))
title('曲率变化率')