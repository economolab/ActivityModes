function plotAllModesEarlyMove(rez,removeEarly,ev,alignEv,plt)
% Function for plotting all modes calculated with and without early
% movement trials

% get field names for each mode
[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');

psth = rez.psth - mean(rez.psth,1);

sample = mode(ev.sample) - mode(ev.(alignEv));
delay  = mode(ev.delay) - mode(ev.(alignEv));
zeroEv  = 0; % corresponds to either gocue or moveonset

% plot each mode
fig = figure;
for i = 1:numel(fns)
    subplot(2,1,i);
    hold on
    for j = 1:numel(plt.conditions)
        cond = plt.conditions(j);
        latent = mySmooth(psth(:,:,cond)*rez.(fns{i}),plt.smooth);
        plot(rez.time,latent,'Color',plt.colors{j},'LineWidth',plt.lw(j),'LineStyle',plt.ls(1));
    end
    hold on;
    for j = 1:numel(plt.conditions)
        cond = plt.conditions(j);
        latent = mySmooth(psth(:,:,cond)*removeEarly.(fns{i}),plt.smooth);
        plot(rez.time,latent,'Color',plt.colors{j+2},'LineWidth',plt.lw(j),'LineStyle',plt.ls(2));
    end
    
    xline(0,'k--','LineWidth',0.5);
    title(fns{i}(1:end-5))
    xlim([-2.5,2.1])
    hold off
end

leg = legend(plt.legend);
leg.Position = [0.60,0.066,0.28,0.13];
sgtitle(plt.title)

% h=axes(fig,'visible','off');
% h.XLabel.Visible='on';
% h.YLabel.Visible='on';
% ylabel(h,'Activity (a.u.)');
% 
% xlabel(h,['Time (s) from ' alignEv]);

end % plotAllModes





