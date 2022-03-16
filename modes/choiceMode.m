function choicemode = choiceMode(obj,meta,cond,epoch,alignEvent,RemoveEarly)
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

 % If early movement trials were identified, exclude them from the
 % trials used to find the mode
 if strcmp(RemoveEarly,'yes')
     trials.ix(obj.earlyMoveix,:) = 0;
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

epochMean = getEpochMean(obj,epochix,trials,meta,RemoveEarly);

[mu,sd] = getEpochStats(epochMean,meta,trials);

% calculate mode according to definition
choicemode = (mu(:, 1) - mu(:, 2))./sqrt(sum(sd(:, 1:2).^2, 2));
% choicemode = ((mu(:,1)-mu(:,3)) + (mu(:,4)-mu(:,2)))./ sqrt(sum(sd.^2,2));
choicemode(isnan(choicemode)) = 0;
choicemode = choicemode./sum(abs(choicemode)); % (ncells,1)

end % choiceMode