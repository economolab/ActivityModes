function contextmode = contextModeMulti(objs,meta,cond,epoch,alignEvent)
% 8. context mode: can identify multiple for different epochs
%      ((hit2AFC - hitAW) / sqrt(sum(sd for each tt ^2));

mu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    % which trials to use for each condition used for finding the mode
    obj = objs{i};
    trials = getTrialsForModeID(objs{i},cond);
    if isfield(obj,'earlyMoveix')
        trials.ix(obj.earlyMoveix,:) = 0;
    elseif isfield(obj,'earlyMoveTrial')
        trials.ix(obj.earlyMoveTrial,:) = 0;
    end
    
    % find time in each trial corresponding to epoch
    epochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        epochix(trix,:) = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch,trix,alignEvent); % (idx1,idx2)
        %epochix(trix,:) = findedges(objs{i}.time,objs{i}.bp,meta(i).dt,epoch,trix,alignEvent); % (idx1,idx2)
    end
    
    epochMean = getEpochMean(objs{i},epochix,trials,meta(i));
    
    [mu{i},sd{i}] = getEpochStats(epochMean,meta(i),trials);
end

mu = cell2mat(mu);
sd = cell2mat(sd);

% calculate mode according to definition
contextmode = (mu(:,1)-mu(:,2))./ sqrt(sum(sd.^2,2));
contextmode(isnan(contextmode)) = 0;
contextmode(isinf(contextmode)) = 0;
contextmode = contextmode./sum(abs(contextmode)); % (ncells,1)


end % actionMode