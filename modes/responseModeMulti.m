function responsemode = responseModeMulti(objs,meta,cond,epoch,alignEvent,rezpsth,psthcond,RemoveEarly)
% response mode: 
%    a. find eigenvectors of baseline subtracted PSTHs using SVD
%       aa. matrix was of size (n x (2t)), where left and right trials concatenated
%       in time
%    b. response mode = eigenvector with largest eigenvalue


mu = cell(numel(objs),1);
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
    
    [mu{i},~] = getEpochStats(epochMean,meta(i),trials);
end

mu = cell2mat(mu);

% get psth, subtract off baseline, concatenate
psth = rezpsth(:,:,psthcond);
for i = 1:numel(psthcond)
    psth(:,:,i) = (squeeze(psth(:,:,i))' - mu(:,i))';
end

X = [psth(:,:,1) ; psth(:,:,2)];

% SVD
% [U,S,V] = svd(X-mean(X)); 
V = myPCA(X - mean(X));
responsemode = V(:,1); % S returned in decreasing order

end % responseMode