function rampingmode = rampingMode(obj,meta,cond,epoch,alignEvent,RemoveEarly)
% ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

 % If early movement trials were identified, exclude them from the
 % trials used to find the mode
if strcmp(RemoveEarly,'yes')
        trials.ix(obj.earlyMoveix,:) = 0;
end

% find time in each trial corresponding to epoch
sampepochix = nan(obj.bp.Ntrials,2);
delayepochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    if strcmp(alignEvent,'firstLick')                   % If everything is aligned to the first lick, use the find edges function that accounts for that
        sampepochix(trix,:)  = findedges_FirstLick(obj.time,obj.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
        delayepochix(trix,:) = findedges_FirstLick(obj.time,obj.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
    else
        sampepochix(trix,:)  = findedges(obj.time,obj.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
        delayepochix(trix,:) = findedges(obj.time,obj.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
    end
end

sampEpochMean  = getEpochMean(obj,sampepochix,trials,meta);
delayEpochMean = getEpochMean(obj,delayepochix,trials,meta);

[sampmu,sampsd] = getEpochStats(sampEpochMean,meta,trials);
[delaymu,delaysd] = getEpochStats(delayEpochMean,meta,trials);
sd =  [sampsd delaysd]; % (clu,cond)

rampingmode = (sampmu - delaymu) ./ sqrt(sum(sd.^2,2));
rampingmode(isnan(rampingmode)) = 0;
rampingmode = rampingmode./sum(abs(rampingmode)); % (ncells,1)


end % outcomeMode