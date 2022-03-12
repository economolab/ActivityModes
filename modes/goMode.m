function gomode = goMode(obj,meta,cond,epoch,alignEvent,RemoveEarly)
% go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% If early movement trials were identified, exclude them from the
 % trials used to find the mode
 if strcmp(RemoveEarly,'yes')
     trials.ix(obj.earlyMoveix,:) = 0;
 end


% find time in each trial corresponding to epoch
postepochix = nan(obj.bp.Ntrials,2);
preepochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    if strcmp(alignEvent,'firstLick')
        postepochix(trix,:)  = findedges_firstLick(obj.time,obj.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
        preepochix(trix,:) = findedges_firstLick(obj.time,obj.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
    else
        postepochix(trix,:)  = findedges(obj.time,obj.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
        preepochix(trix,:) = findedges(obj.time,obj.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
    end
end

postEpochMean  = getEpochMean(obj,postepochix,trials,meta,RemoveEarly);
preEpochMean = getEpochMean(obj,preepochix,trials,meta,RemoveEarly);

[postmu,postsd] = getEpochStats(postEpochMean,meta,trials);
[premu,presd] = getEpochStats(preEpochMean,meta,trials);
sd =  [postsd presd]; % (clu,cond)

gomode = (postmu - premu) ./ sqrt(sum(sd.^2,2));
gomode(isnan(gomode)) = 0;
gomode = gomode./sum(abs(gomode)); % (ncells,1)


end % outcomeMode