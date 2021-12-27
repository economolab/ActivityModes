function obj = getSeq(obj,meta)

edges = meta.tmin:meta.dt:meta.tmax;                % Times -3 s before go-cue and 3s after go-cue
obj.time = edges + meta.dt/2;                       % Shift edges by half a time bin
obj.time = obj.time(1:end-1);

% get psths by condition
obj.psth = zeros(numel(obj.time),numel(meta.cluid),numel(obj.condition));
for i = 1:numel(meta.cluid)                         % For each cluster...
    curClu = meta.cluid(i);                         
    for j = 1:numel(obj.condition)                  % For each behavioral condition...
        trix = meta.trialid{j};                     % Get the trial indices that correspond to the jth behavioral condition 
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);   % Returns a logical array of length(#spikes for ith cluster)
                                                                     % where true values are spikes that occur in a trial of the jth condition   
        
        % Obj.clu{meta.probe}(curClu).trialtm_aligned(spkix) = The aligned
        % spike times for ith cluster within each jth trial
        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);    % Finds the number of spikes for this cluster that occur within each time bin (histogram)
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
obj.trialpsth = zeros(numel(obj.time),numel(meta.cluid),obj.bp.Ntrials);
spiketrains = zeros(numel(obj.time),numel(meta.cluid),obj.bp.Ntrials);      % For MA881 project
for i = 1:numel(meta.cluid)
    curClu = meta.cluid(i);
    for j = 1:obj.bp.Ntrials
                
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, j);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        spiketrains(:,i,j) = N;                                             % For MA881 project
                                                                            % # of time bins x # of neurons # of trials
        obj.trialpsth(:,i,j) = mySmooth(N./meta.dt,meta.smooth);

    end
end

% For MA881 project
cond_spiketrain = cell(1,7);
cond_spiketrain{1,7} = [];
for j = 1:numel(obj.condition)                  % For each behavioral condition...
        trix = meta.trialid{j};  
        cond_spiketrain{j} = spiketrains(:,:,trix);
end 

end % getSeq