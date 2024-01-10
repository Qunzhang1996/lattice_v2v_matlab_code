function [allTS, allDT, numTS] = SamplingEndcontions(refPath, laneWidth, laneBounds, egoState, curActorState, safetyGap, speedLimit, timeHorizons)
    % Generate cruise control states.
    [termStatesCC,timesCC] = SamplingBasicCruiseControl(...
        refPath,laneWidth,laneBounds,egoState,speedLimit,timeHorizons);
    
    % Generate lane change states.
    [termStatesLC,timesLC] = SamplingBasicLaneChange(...
        refPath,laneWidth,laneBounds,egoState,timeHorizons);

    % Generate vehicle following states.
    [termStatesF,timesF] = SamplingBasicLeadVehicleFollow(...
        refPath,laneWidth,laneBounds,safetyGap,egoState,curActorState,timeHorizons);
    
    % Combine the terminal states and times.
    allTS = [termStatesCC; termStatesLC; termStatesF];
    allDT = [timesCC; timesLC; timesF];
    numTS = [numel(timesCC); numel(timesLC); numel(timesF)];
end