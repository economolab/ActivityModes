function choicemode = choiceMode(obj,meta,cond,epoch,alignEvent,RemoveShort)
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

 % Remove trials with delay lengths less than 0.6 from identifying the
 % choice mode
 shortix = [];
 for c = 1:numel(cond)
     if isfield(meta,'del_trialid')
        shortix = [shortix;find(meta.del_trialid{c}==1 | meta.del_trialid{c}==2)];
     end
 end
 if strcmp(RemoveShort,'yes')
     trials.ix(shortix,:) = 0;
 end

 % find time in each trial corresponding to epoch
 epochix = nan(obj.bp.Ntrials,2);
 for trix = 1:obj.bp.Ntrials
     if strcmp(alignEvent,'firstLick')                   % If everything is aligned to the first lick, use the find edges function that accounts for that
         epochix(trix,:) = findedges_FirstLick(obj.time,obj.bp,epoch,trix,alignEvent); % (idx1,idx2)
     else                                                % Otherwise, use the normal find edges function
         epochix(trix,:) = findedges(obj.time,obj.bp,epoch,trix,alignEvent); % (idx1,idx2)
     end
 end

epochMean = getEpochMean(obj,epochix,trials,meta,shortix,RemoveShort);

[mu,sd] = getEpochStats(epochMean,meta,trials);

% calculate mode according to definition
choicemode = (mu(:, 1) - mu(:, 2))./sqrt(sum(sd(:, 1:2).^2, 2));            % Delay period coding dimension
% choicemode = ((mu(:,1)-mu(:,3)) + (mu(:,4)-mu(:,2)))./
% sqrt(sum(sd.^2,2));                                                       % Nuo Li's definition
choicemode(isnan(choicemode)) = 0;
choicemode = choicemode./sum(abs(choicemode)); % (ncells,1)

end % choiceMode