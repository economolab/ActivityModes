function stimmode = stimulusModeMulti(objs,meta,cond,epoch,alignEvent,RemoveEarly)
% 1. stimulus mode: defined during stimulus (sample) period
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));

mu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    % which trials to use for each condition used for finding the mode
    trials = getTrialsForModeID(objs{i},cond);
    
    % If early movement trials were identified, exclude them from the
    % trials used to find the mode
    if strcmp(RemoveEarly,'yes')
     trials.ix(objs{i}.earlyMoveix,:) = 0;
    end
    
    % find time in each trial corresponding to epoch
    epochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        if strcmp(alignEvent,'firstLick')                   % If everything is aligned to the first lick, use the find edges function that accounts for that
            epochix(trix,:) = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch,trix,alignEvent); % (idx1,idx2)
        else                                                % Otherwise, use the normal find edges function
            epochix(trix,:) = findedges(objs{i}.time,objs{i}.bp,meta(i).dt,epoch,trix,alignEvent); % (idx1,idx2)
        end
    end
    
    epochMean = getEpochMean(objs{i},epochix,trials,meta(i),RemoveEarly);
    
    [mu{i},sd{i}] = getEpochStats(epochMean,meta(i),trials);
end

mu = cell2mat(mu);
sd = cell2mat(sd);

% calculate mode according to definition
stimmode = ((mu(:,1)-mu(:,4)) + (mu(:,3)-mu(:,2)))./ sqrt(sum(sd.^2,2));
stimmode(isnan(stimmode)) = 0;
stimmode = stimmode./sum(abs(stimmode)); % (ncells,1)

end % stimulusMode





