function actionmode = actionMode(obj,meta,cond,epoch,alignEvent,RemoveEarly)
% action mode: defined during mvmt init (0.1 to 0.3 s rel to go cue)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% % get move onset times
% if ~isfield(obj.bp.ev,'moveOnset') && strcmpi(epoch,'moveonset')
%     obj = findMoveOnset(obj);
% end

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

% calculate mode according to definition
actionmode = (mu(:,1)-mu(:,2))./ sqrt(sum(sd.^2,2));
actionmode(isnan(actionmode)) = 0;
actionmode = actionmode./sum(abs(actionmode)); % (ncells,1)


end % actionMode