clear,clc,close all

if ispc
    pth = 'C:\Code\activity_modes';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/activity_modes';
end

addpath(genpath(pwd))

% find all 7 activity modes as described in :
% Thalamus drives diverse responses in frontal cortex during decision-making
% Weiguo Yang, Sri Laasya Tipparaju, Guang Chen, Nuo Li

% 1. stimulus mode: defined during stimulus (sample) period
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));
% 3. action mode: defined during mvmt init (0.1 to 0.3 s rel to go cue)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));
% 4. outcome mode: defined during response epoch (0 to 1.3 s rel go cue)
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));
% 5. ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));
% 6. go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));
% 7. response mode: 
%    a. find eigenvectors of basline subtracted PSTHs using SVD
%       aa. matrix was of size (n x (2t)), where left and right trials concatenated
%       in time
%    b. response mode = eigenvector with largest eigenvalue


%% TODO
% subtract out modes found using 2afc from aw context. See what's left.
% preprocess data other than normalize???
% handle multiple probes

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'goCue' or 'moveOnset'

params.lowFR               = 1; % remove clusters firing less than this val

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&~autowater&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&~autowater&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&~autowater&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&~autowater&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&~autowater&~early'}; % ignore
params.condition(8) = {'hit&~stim.enable&~autowater&~early'}; % hit 2afc
params.condition(9) = {'hit&~stim.enable&autowater&~early'}; % hit aw

% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};


params.lowFR               = 1;      % when doing single trial analysis, remove clusters with avg firing rate < params.lowFR

params.probe = 1;
params.probeArea = 'ALM';

params.tmin = -2.5;
params.tmax = 1.5;
params.dt = 1/200;

params.smooth = 15;

params.quality = {'all'};

%% SET METADATA
% experiment meta data
meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.anm = 'JEB7';
meta.date = '2021-04-29';
% meta.datafn = findDataFn(meta);
meta.datafn = 'data_structure_JEB7_2021-04-29.mat';

%% LOAD DATA

[meta,params,obj] = loadAndProcessData(meta,params);

%% label move or non-move
[obj.movix,obj.movtime] = getMoveIdx(obj,params);

%% ACTIVITY MODES
rez.time = obj.time;
rez.psth = obj.psth;
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;

%% context mode
% hit2afc - hitaw during presample period
cond{1} = params.condition{8}; % hit 2afc
cond{2} = params.condition{9}; % hit aw
rez.context_mode = (obj.presampleFR(:,8) - obj.presampleFR(:,9)) ./ sqrt(sum(obj.presampleSigma(:,8:9).^2,2));
clear cond

%% stimulus mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'sample';
rez.stimulus_mode = stimulusMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% choice mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode = choiceMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% action mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'action';
rez.action_mode = actionMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% outcome mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'outcome';
rez.outcome_mode = outcomeMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% ramping mode
cond{1} = params.modecondition{5};
epoch = {'presample','delay'};
rez.ramping_mode = rampingMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% go mode
cond{1} = params.modecondition{5};
epoch = {'postgo','prego'};
rez.go_mode = goMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% response mode 
% cond{1} = params.modecondition{1};
% cond{2} = params.modecondition{2};
% psthcond = [1,2];
% epoch = 'presample'; % used to estimate baseline firing rate
% rez.response_mode = responseMode(obj,params,cond,epoch,rez.alignEvent,psthcond);
% clear cond


%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(params.cluid),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end

%% remainder mode
modesToKeep = eye(size(obj.psth,2)) - (orthModes*orthModes');

residualpsth = nan(size(obj.psth));
for i = 1:size(obj.psth,3)
    residualpsth(:,:,i) = obj.psth(:,:,i) * modesToKeep;
end

X = [residualpsth(:,:,1) ; residualpsth(:,:,2)]; % left and right 2afc

% SVD
V = myPCA(X - mean(X));
rez.remainder1_mode = V(:,1); % S returned in decreasing order
rez.remainder2_mode = V(:,2); % S returned in decreasing order
rez.remainder3_mode = V(:,3); % S returned in decreasing order


%% PLOTS

% MODES VIZ

% % plot correct trials alone
% plt.title = 'Correct Trials';
% plt.legend = {'Right Hit','Left Hit'};
% plt.conditions = [1,2];
% plt.lw = [2 2];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0]};
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 
% 
% % plot correct trials and error trials
% plt.title = 'Correct and Error Trials';
% plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
% plt.conditions = [1,2,3,4];
% plt.lw = [2.5 2.5 1.5 1.5];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0], ...
%                  [0.5 0.5 1],[1 0.5 0.5]};
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% plot correct trials and AW trials
plt.title = '2AFC and Autowater (Correct) Trials';
plt.legend = {'Right 2AFC','Left 2AFC','Right AW', 'Left AW'};
plt.conditions = [1,2,5,6];
plt.lw = [2.7 2.7 2.7 2.7];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0], ...
                 [190, 3, 252]./255,[252, 190, 3]./255};
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% % plot correct and ignore trials
% plt.title = 'Correct and Ignore Trials';
% plt.legend = {'Right Hit','Left Hit','Ignore'};
% plt.conditions = [1,2,7];
% plt.lw = [2 2 2];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0],[0 0 0]};
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% figure; 
% psth_jaw = nanmean(obj.trialpsth(:,:,movtrials),3);
% latent = mySmooth(psth_jaw*rez.jaw_mode,30);
% plot(obj.time,latent,'Color','g','LineWidth',2); hold on
% psth_nojaw = nanmean(obj.trialpsth(:,:,nonmovtrials),3);
% latent = mySmooth(psth_nojaw*rez.jaw_mode,30);
% plot(obj.time,latent,'Color','r','LineWidth',2); hold off
% title('Jaw Mode')
% legend('movement trials','non-movement trials')
% ax = gca;
% ax.FontSize = 20;




% ORTHOGONALITY VIZ
% dotProductModes(rez,modes,'NOT ORTHOGONALIZED')
% dotProductModes(rez,orthModes,'ORTHOGONALIZED')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
fn = fn{1};

end % loadRawDataObj












