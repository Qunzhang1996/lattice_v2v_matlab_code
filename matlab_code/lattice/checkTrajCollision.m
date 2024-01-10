function isColliding = checkTrajCollision(egoPoses, actorPoses, Geometry)
isColliding = true;
for i = 1:length(actorPoses)
    if checkOneTrajCollision(egoPoses, actorPoses(i), Geometry)
        return;
    end
end
isColliding = false;
end

function isColliding = checkOneTrajCollision(egoPoses, actorPoses, Geometry)
% Geometry.Length = carLen; % in meters
% Geometry.Radius = carWidth/2; % in meters
% Geometry.FixedTransform = -carRear; % in meters
len1 = size(egoPoses.States,1);
len2 = size(actorPoses.States,1);
p1 = zeros(2, 31);
v1 = zeros(2, 31);
p2 = zeros(2, 31);
v2 = zeros(2, 31);
D1 = Geometry.Length;
R1 = Geometry.Radius;
D2 = D1;
R2 = R1;
for i = 1:31
    if i > len1
        v1(:,i) = v1(:,i-1);
        p1(:,i) = p1(:,i-1);
    else
        v1(:,i) = [cos(egoPoses.States(i,3)); sin(egoPoses.States(i,3))];
        p1(:,i) = egoPoses.States(i,1:2)'+Geometry.FixedTransform*v1(:,i);
    end
end
for i = 1:31
    if i > len2
        v2(:,i) = v2(:,i-1);
        p2(:,i) = p2(:,i-1);
    else
        v2(:,i) = [cos(actorPoses.States(i,3)); sin(actorPoses.States(i,3))];
        p2(:,i) = actorPoses.States(i,1:2)'+Geometry.FixedTransform*v2(:,i);
    end
end
[isCollidings, ~] = checkCollisionCapsule(p1, v1, D1, R1, p2, v2, D2, R2);
isColliding = any(isCollidings);
end
