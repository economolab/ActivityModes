function obj = alignSpikes(obj,meta,params)

if strcmpi(params.alignEvent,'moveOnset')
    obj = findMoveOnset(obj); % assigns obj.bp.ev.moveOnset
end

if strcmp(params.alignEvent,'firstLick')
    % get first lick time for left and right licks
    temp = obj.bp.ev.lickL;
    idx = ~cellfun('isempty',temp);
    outL = zeros(size(temp));
    outL(idx) = cellfun(@(v)v(1),temp(idx));
    temp = obj.bp.ev.lickR;
    idx = ~cellfun('isempty',temp);
    outR = zeros(size(temp));
    outR(idx) = cellfun(@(v)v(1),temp(idx));
    firstLick = zeros(size(temp));
    % firstLick = min(outL,outR), except when outL||outR == 0
    outL(outL==0) = nan;
    outR(outR==0) = nan;
    firstLick = nanmin(outL,outR);
    firstLick(isnan(firstLick)) = 0;
    obj.bp.ev.(params.alignEvent) = firstLick;
end

% align spikes to params.alignEvent
for clu = 1:numel(obj.clu{meta.probe})
    tmp = obj.clu{meta.probe}(clu).trial;                        % To account for there being more SpikeGLX trials than Bpod trials...                     
    tmp(tmp > obj.bp.Ntrials) = [];                              % Remove spikes that occurred in a SpikeGLX trial that was not registered by Bpod
    event = obj.bp.ev.(params.alignEvent)(tmp);                  % Finds the timing of the aligning event for each trial that a spike from this cluster occurs in
    trialtm = obj.clu{meta.probe}(clu).trialtm;
    trialtm = trialtm(1:length(event));                          % Remove spikes that occurred in a SpikeGLX trial that was not registered by Bpod 
    obj.clu{meta.probe}(clu).trialtm_aligned = trialtm - event;
end

end % alignSpikes