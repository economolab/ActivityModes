function trials = getTrialsForModeID(obj,cond)
varnames = getStructVarNames(obj);
for i = 1:numel(varnames)
    eval([varnames{i} ' = obj.bp.' varnames{i} ';']);
    
    if eval(['numel(' varnames{i} ')==obj.bp.Ntrials && isrow(' varnames{i} ')'])
        eval([varnames{i} '=' varnames{i} ''';']);
    end
end

trials.N = numel(cond);
trials.ix = zeros(obj.bp.Ntrials, trials.N);
for i = 1:numel(cond)
    curfilt = cond{i};
    trials.ix(:,i) = eval(curfilt);
end
end % getTrials
