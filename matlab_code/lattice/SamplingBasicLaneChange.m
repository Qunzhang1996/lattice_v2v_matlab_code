function [terminalStates, times] = SamplingBasicLaneChange(refPath, laneWidth,laneBounds, egoState, dt)   
    if egoState(5) == 0
        terminalStates = [];
        times = [];
    else
        % Convert ego state to Frenet coordinates
        frenetState = global2frenet(refPath, egoState);

        % Get current lane
        curLane = PredictLane(frenetState, laneWidth,laneBounds, 0);

        % Determine if future lanes are available
        adjacentLanes = curLane+[-1 1];
        validLanes = adjacentLanes > 0 & adjacentLanes <= 4;

        % Calculate lateral deviation for adjacent lanes
        lateralOffset = (2-adjacentLanes(validLanes)+.5)*laneWidth;
        numLane = nnz(validLanes);

        % Calculate terminal states
        terminalStates = zeros(numLane*numel(dt),6);
        terminalStates(:,1) = nan;
        terminalStates(:,2) = egoState(5);
        terminalStates(:,4) = repelem(lateralOffset(:),numel(dt),1);
        times = repmat(dt(:),numLane,1);
    end
end