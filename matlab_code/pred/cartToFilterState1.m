function filterState = cartToFilterState1(cartState, refPath)
% Assemble as global state to use global2frenet
numStates = size(cartState,2);
x = cartState(1,:);
y = cartState(3,:);
theta = atan2(cartState(4,:),cartState(2,:));
speed = sqrt(cartState(4,:).^2 + cartState(2,:).^2);
theta(speed <= 2) = 0;
kappa = zeros(1,numStates);
accel = zeros(1,numStates);
globalState = [x;y;theta;kappa;speed;accel];
% Convert to Frenet state
% refPath = helperGetReferencePath;
try
    frenetState = global2frenet(refPath,globalState')';
catch
    globalState(3) = globalState(3) + pi;
    frenetState = global2frenet(refPath,globalState')';
end

% Convert to filter state
s = frenetState(1,:);
ds = frenetState(2,:);
d = frenetState(4,:);
dd = frenetState(5,:).*ds;
filterState = [s;ds;d;dd];

end