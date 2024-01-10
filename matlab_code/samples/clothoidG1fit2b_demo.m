% clothoidG1fit2b_demo.m

waypoint1s = [0 0 -0.3531; ...
            50 -10 0.1148; ...
            100 10 0.4193; ...
            150 20 -0.0685
            200 10 -0.1484
            250 10 0.0742];
waypoint2s = [0 0 0; ...
            50 -10 0; ...
            100 10 0; ...
            150 20 0
            200 10 0
            250 10 0];
waypoints = waypoint1s;
n = size(waypoints,1);
initialPosition = complex(waypoints(:,1), waypoints(:,2));
course = waypoints(:,3);
[initialCurvature, finalCurvature, arcLengths] = clothoidG1fit2(initialPosition(1:n-1),course(1:n-1),initialPosition(2:n),course(2:n));
cumulativeLength = [0;cumsum(arcLengths)];
segStarts1 = [real(initialPosition) imag(initialPosition) robotics.internal.wrapToPi(course)];
segStarts2 = [[initialCurvature; finalCurvature(end)], ...
              [(finalCurvature-initialCurvature)./arcLengths; (finalCurvature(end)-initialCurvature(end))/arcLengths(end)],...
              cumulativeLength];
segStarts  = [segStarts1, segStarts2];
hold on;
ss  = [];
ths = [];
ks  = [];
for i = 1:(size(segStarts,1)-1)
    plot(segStarts(i,1),segStarts(i,2), 'o');
    segs = 0:0.1:(segStarts(i+1,6)-segStarts(i,6));
    [x,y,th,k] = FrenetReferencePath.clothoid(segStarts(i,1),...
                                              segStarts(i,2),...
                                              segStarts(i,3),...
                                              segStarts(i,4),...
                                              segStarts(i,5),...
                                              segs);
    ss  = [ss,segs+segStarts(i,6)];
    ths = [ths,th];
    ks  = [ks,k];
    plot(x,y, 'marker', '.');
end
figure
subplot(2,1,1)
plot(ss, ths)
title('航向角')
subplot(2,1,2)
plot(ss, ks)
title('曲率')
