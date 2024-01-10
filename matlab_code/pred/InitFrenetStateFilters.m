function filter = InitFrenetStateFilters(cartState, refPath)
% Convert the state to Frenet coordinates
state = cartToFilterState(cartState, refPath); 
% Define state covariance
stateCov = blkdiag(5,100,4,5);
% Define filter
filter = trackingEKF(@stateTransFcn,@measurementFcn,state,...
    'StateTransitionJacobianFcn',@stateTransJacobianFcn,...
    'StateCovariance',stateCov,...
    'MeasurementNoise',blkdiag(1,1,1,1),...
    'HasAdditiveProcessNoise',false,...
    'ProcessNoise',blkdiag(1,100));
end

function state = stateTransFcn(state, w, dT)
if isscalar(w)
    w = repmat(w,[4 size(state,2)]);
end
tau = 3;
state(1,:) = state(1,:) + state(2,:)*dT + w(1,:)*dT^2/2;
state(2,:) = state(2,:) + w(1,:)*dT;
state(3,:) = state(3,:) + state(4,:)*tau*(1 - exp(-dT/tau)) + w(2,:)*dT^2/2;
state(4,:) = state(4,:)*exp(-dT/tau) + w(2,:)*dT;
end

function [dfdx, dfdw] = stateTransJacobianFcn(~, ~, dT)
    tau = 3;
    dfdx = [1 dT 0 0;0 1 0 0;0 0 1 tau*(1 - exp(-dT/tau));0 0 0 exp(-dT/tau)];
    dfdw = [dT^2/2 0;dT 0;0 dT^2/2;0 dT];
end

function z = measurementFcn(x, refPath)
    z = filterToCartState(x, refPath);
end

function filterState = cartToFilterState(cartState, refPath)
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

function cartState = filterToCartState(filterState, refPath)

% Assemble as Frenet state to use frenet2global
numStates = size(filterState,2);
s = filterState(1,:);
ds = filterState(2,:);
dds = zeros(1,numStates);
d = filterState(3,:);
dd = filterState(4,:);
ddbyds = zeros(1,numStates);
ddbyds(abs(ds) > 0.5) = dd(abs(ds) > 0.5)./ds(abs(ds) > 0.5);
dd2ds2 = zeros(1,numStates);
frenetState = [s;ds;dds;d;ddbyds;dd2ds2];

% Convert to global state
% refPath = helperGetReferencePath;
globalState = frenet2global(refPath,frenetState')';

% Convert to cartesian state
x = globalState(1,:);
y = globalState(2,:);
speed = globalState(5,:);
theta = globalState(3,:);
vx = speed.*cos(theta);
vy = speed.*sin(theta);
cartState = [x;vx;y;vy];

end