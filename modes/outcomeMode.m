function outcomemode = outcomeMode(obj,meta,cond,epoch,alignEvent,RemoveEarly)
% outcome mode: defined during response epoch (0 to 1.3 s rel go cue)
%       ((hitR - missR) + (hitL - missL)) / sqrt(sum(sd for each tt ^2));

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

epochMean = getEpochMean(obj,epochix,trials,meta);

[mu,sd] = getEpochStats(epochMean,meta,trials);

outcomemode = ((mu(:,1)-mu(:,3)) + (mu(:,2)-mu(:,4)))./ sqrt(sum(sd.^2,2));
outcomemode(isnan(outcomemode)) = 0;
outcomemode = outcomemode./sum(abs(outcomemode)); % (ncells,1)


end % outcomeMode