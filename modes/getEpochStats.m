function [mu,sd] = getEpochStats(epochMean,meta,trials)
% calculate mu and stdev across trials
mu = nan(numel(meta.cluid),trials.N);
sd = nan(size(mu));
for cluix = 1:numel(meta.cluid) % for each cluster
    for cnd = 1:trials.N
        mu(cluix,cnd) = nanmean(epochMean(cluix,:,cnd));
        sd(cluix,cnd) = nanstd(epochMean(cluix,:,cnd));
    end
end
end % getEpochStats