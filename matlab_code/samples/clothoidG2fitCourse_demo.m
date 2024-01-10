% clothoidG2fitCourse_demo.m

waypoints = [0 0 -0.3531; ...
            50 -10 0.1148; ...
            100 10 0.4193; ...
            150 20 -0.0685
            200 10 -0.1484
            250 10 0.0742];
course = clothoidG2fitCourse(waypoints(:,1:2));
disp(course)