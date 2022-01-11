function trialNums = findTrials(obj, conditions)

if ~isfield(obj.bp.autowater, 'nums')           %If 'nums' is not a field in obj.bp.autowater...(aka behavior run on MasterProtocol_v6)
    tmp = obj.bp.autowater;
    obj.bp = rmfield(obj.bp, 'autowater');      %Remove the 'autowater' field from the object
    obj.bp.autowater.nums = tmp + (tmp-1)*-2;   %And re-create it such that it contains autowater.nums==1 and 2 (as in old MasterProtocol_v4)
end

varnames = getStructVarNames(obj);
for i = 1:numel(varnames)
    eval([varnames{i} ' = obj.bp.' varnames{i} ';']);
    
    if eval(['numel(' varnames{i} ')==obj.bp.Ntrials && isrow(' varnames{i} ')'])
        eval([varnames{i} '=' varnames{i} ''';']);
    end
end

mask = zeros(obj.bp.Ntrials, numel(conditions));
for i = 1:numel(conditions)
    mask(:,i) = eval(conditions{i}); 
    trialNums{i} = find(mask(:,i));
end

end % findTrials