function contextmode = contextMode(obj,meta,cond,epoch,alignEvent)
% 2. context mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,meta.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

epochMean = getEpochMean(obj,epochix,trials,meta);

[mu,sd] = getEpochStats(epochMean,meta,trials);

% calculate mode according to definition
contextmode = ((mu(:,1)-mu(:,3)) + (mu(:,4)-mu(:,2)))./ sqrt(sum(sd.^2,2));
contextmode(isnan(contextmode)) = 0;
contextmode = contextmode./sum(abs(contextmode)); % (ncells,1)

end % contextMode