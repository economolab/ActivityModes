function obj = getSeq(obj,meta)

edges = meta.tmin:meta.dt:meta.tmax;
obj.time = edges + meta.dt/2;
obj.time = obj.time(1:end-1);

% get psths by condition
obj.psth = zeros(numel(obj.time),numel(meta.cluid),numel(obj.condition));
for i = 1:numel(meta.cluid)
    curClu = meta.cluid(i);
    for j = 1:numel(obj.condition)
        trix = meta.trialid{j};
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
obj.trialpsth = zeros(numel(obj.time),numel(meta.cluid),obj.bp.Ntrials);
for i = 1:numel(meta.cluid)
    curClu = meta.cluid(i);
    for j = 1:obj.bp.Ntrials
                
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, j);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialpsth(:,i,j) = mySmooth(N./meta.dt,meta.smooth);

    end
end


end % getSeq