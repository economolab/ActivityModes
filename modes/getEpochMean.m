function epochMean = getEpochMean(obj,epochix,trials,meta,shortix,RemoveShort)
% slice data by trials and epochix
% compute the mean during the epoch for each trial, cluster
nTrials = max(sum(trials.ix));
nTime = max(diff(epochix'));
psth = nan(nTime,numel(meta.cluid),nTrials,trials.N); % (time,clu,trials,cond)
epochMean = nan(numel(meta.cluid),nTrials,trials.N);  % (cells x trials x cond needed for mode)
for cluix = 1:numel(meta.cluid) % for each cluster
    for cnd = 1:trials.N        % for each cond
        cndtrid = find(trials.ix(:,cnd));   
        for trix = 1:numel(cndtrid)     % for every trial in cnd
            trid = cndtrid(trix);
            
            if strcmp(RemoveShort,'yes')        % If you are removing early trials from the mode calculation...
                if ~ismember(trid,shortix)      % Only do the following if the current trial is not an early move trial
                    e1 = epochix(trid,1);               % Get the index that corresponds to the beginning of the given epoch in curr trial
                    e2 = epochix(trid,2);               % Get the index that corresponds to the end of the given epoch in curr trial
                    
                    psth(1:diff(epochix(trid,:))+1,cluix,trix,cnd) = ...    % For the current trial, get the cell's PSTH during the specified epoch 
                        obj.trialpsth(e1:e2,cluix,trid);
                    
                    % calculate the avg firing rate during the epoch for trial
                    epochMean(cluix,trix,cnd) = mean(psth(:,cluix,trix,cnd),1,'omitnan');   % For the current trial, get the cell's avg FR during the specified epoch
                end
           
            else
                e1 = epochix(trid,1);
                e2 = epochix(trid,2);
                
                psth(1:diff(epochix(trid,:))+1,cluix,trix,cnd) = ...
                    obj.trialpsth(e1:e2,cluix,trid);
                
                % calculate the avg firing rate during the epoch for trial
                epochMean(cluix,trix,cnd) = mean(psth(:,cluix,trix,cnd),1,'omitnan');
            end
        end
    end
end

end % getEpochMean