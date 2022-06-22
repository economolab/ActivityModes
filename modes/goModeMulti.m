function gomode = goModeMulti(objs,meta,cond,epoch,alignEvent,RemoveShort)
% go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));

postmu = cell(numel(objs),1);
premu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    obj = objs{i};
    met = meta(i);

    % which trials to use for each condition used for finding the mode
    trials = getTrialsForModeID(obj,cond);

    % If early movement trials were identified, exclude them from the
    % trials used to find the mode
    shortix = [];
    for c = 1:numel(cond)
        if isfield(met,'del_trialid')
            shortix = [shortix;find(met.del_trialid{c}==1 | me.del_trialid{c}==2)];
        end
    end
    if strcmp(RemoveShort,'yes')
        trials.ix(shortix,:) = 0;
    end

    % find time in each trial corresponding to epoch
    postepochix = nan(objs{i}.bp.Ntrials,2);
    preepochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        if strcmp(alignEvent,'firstLick')
            postepochix(trix,:)  = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
            preepochix(trix,:) = findedges_FirstLick(objs{i}.time,objs{i}.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
        else
            postepochix(trix,:)  = findedges(objs{i}.time,objs{i}.bp,epoch{1},trix,alignEvent); % (idx1,idx2)
            preepochix(trix,:) = findedges(objs{i}.time,objs{i}.bp,epoch{2},trix,alignEvent); % (idx1,idx2)
        end
    end

    postEpochMean  = getEpochMean(objs{i},postepochix,trials,meta(i),shortix,RemoveShort);
    preEpochMean   = getEpochMean(objs{i},preepochix,trials,meta(i),shortix,RemoveShort);

    [postmu{i},postsd] = getEpochStats(postEpochMean,meta(i),trials);
    [premu{i},presd]   = getEpochStats(preEpochMean,meta(i),trials);
    sd{i} =  [postsd presd]; % (clu,cond)

end

postmu = cell2mat(postmu);
premu = cell2mat(premu);
sd = cell2mat(sd);


gomode = (postmu - premu) ./ sqrt(sum(sd.^2,2));
gomode(isnan(gomode)) = 0;
gomode = gomode./sum(abs(gomode)); % (ncells,1)


end % outcomeMode