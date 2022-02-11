function obj = timeWarp(obj,params)

% get first params.Nlicks lick times, for each trial
[lickStart,lickEnd,lickDur] = findLickTimes(obj,params.nLicks);

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
%find median lick times for each lick across trials
% only using trials where a lick was present in median calculation
med = findMedianLickTimes(lickStart,lickEnd,lickDur, params.nLicks);

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% find fit for each trial and each lick
p = trialWarpFits(lickStart,lickEnd,lickDur,med,obj,params);

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% warp spike times for each cluster (only spike times between go cue and first params.nLicks licks get
% warped)
for cluix = 1:numel(obj.clu{params.probe})
    
    disp(['Warping spikes for cluster ' num2str(cluix) ' / ' num2str(numel(obj.clu{params.probe}))])
    obj.clu{params.probe}(cluix).trialtm_warped = obj.clu{params.probe}(cluix).trialtm;

    for trix = 1:obj.bp.Ntrials
        % find spike times for current trial
        spkmask = ismember(obj.clu{params.probe}(cluix).trial,trix);
        spkix = find(spkmask);
        spktm = obj.clu{params.probe}(cluix).trialtm(spkmask);
        
        
        for lix = 1:numel(lickStart{trix}) % lick index for current trial
            p_cur = p{trix,lix}; % current fit parameters for current trial and lick number
            
            ls = lickStart{trix}(lix); % current trial lick(lix) start times
            le = lickEnd{trix}(lix); % current trial lick(lix) end times
            ld = lickDur{trix}(lix); % current trial lick(lix) durations
            
            % find spike ix b/w ls and le (these are the spikes that will
            % be time warped)
            mask = (spktm>ls) & (spktm<le);
            
            tm = spktm(mask); % spike times between gocue and licks or previous and current lick
            
            % warp
            warptm = polyval(p_cur,tm);
            obj.clu{params.probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
        end
        
    end
    
end


end % timeWarp




























