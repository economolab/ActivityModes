function rampingmode = rampingModeMulti(objs,meta,cond,epoch,alignEvent,RemoveEarly)
% ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));

sampmu = cell(numel(objs),1);
delaymu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    obj = objs{i};
    % which trials to use for each condition used for finding the mode
    trials = getTrialsForModeID(objs{i},cond);
    if strcmp(RemoveEarly,'yes')
        trials.ix(objs{i}.earlyMoveix,:) = 0;
    end
    
    % find time in each trial corresponding to epoch
    sampepochix = nan(objs{i}.bp.Ntrials,2);
    delayepochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        if strcmp(alignEvent,'firstLick')
            sampepochix(trix,:)  = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
            delayepochix(trix,:) = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
        else
            sampepochix(trix,:)  = findedges(objs{i}.time,objs{i}.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
            delayepochix(trix,:) = findedges(objs{i}.time,objs{i}.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
        end

    end
    
    sampEpochMean  = getEpochMean(objs{i},sampepochix,trials,meta(i),RemoveEarly);
    delayEpochMean = getEpochMean(objs{i},delayepochix,trials,meta(i),RemoveEarly);
    
    [sampmu{i},sampsd]   = getEpochStats(sampEpochMean,meta(i),trials);
    [delaymu{i},delaysd] = getEpochStats(delayEpochMean,meta(i),trials);
    sd{i} =  [sampsd delaysd]; % (clu,cond)
end

sampmu = cell2mat(sampmu);
delaymu = cell2mat(delaymu);
sd = cell2mat(sd);

rampingmode = (sampmu - delaymu) ./ sqrt(sum(sd.^2,2));
rampingmode(isnan(rampingmode)) = 0;
rampingmode(isinf(rampingmode)) = 0;
rampingmode = rampingmode./sum(abs(rampingmode)); % (ncells,1)


end % outcomeMode