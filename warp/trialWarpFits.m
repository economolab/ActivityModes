function p = trialWarpFits(lickStart,lickEnd,lickDur,med,obj,params)


p = cell(obj.bp.Ntrials,params.nLicks);
for trix = 1:obj.bp.Ntrials
    ls = lickStart{trix}; % current trial lick start times
    le = lickEnd{trix}; % current trial lick end times
    ld = lickDur{trix}; % current trial lick durations
        
    for lix = 1:numel(ls) % lick index
        % median values for current lick
        mls = med.lickStart(lix); % median lick start time
        mle = med.lickEnd(lix); % median lick end time
        mld = med.lickDur(lix); % median lick duration
        
        % warp lick times
        % - from (lick start time, lick end time) to (median lick start
        % time, median lick end time)
        x = [ls(lix) le(lix)];
        y = [mls mle];

        % fit original time to warped time
        p{trix,lix} = polyfit(x,y,1);
        
    end
    
end

end