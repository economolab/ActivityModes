function plotAllModes(rez,ev,alignEv,plt)

% get field names for each mode
[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');

psth = rez.psth;
datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
eigsum = sum(eig(datacov));

sample = mode(ev.sample) - mode(ev.(alignEv));
delay  = mode(ev.delay) - mode(ev.(alignEv));
zeroEv  = 0; % corresponds to either gocue or moveonset

% plot each mode
fig = figure;
for i = 1:numel(fns)
    subplot(4,3,i);
    hold on
    for j = 1:numel(plt.conditions)
        cond = plt.conditions(j);
        latent = mySmooth(psth(:,:,cond)*rez.(fns{i}),plt.smooth);
        plot(rez.time,latent,'Color',plt.colors{j},'LineWidth',plt.lw(j));
    end
    
    varExp = var_proj(rez.(fns{i}), datacov, eigsum);

    
    xline(sample,'k--','LineWidth',0.5);
    xline(delay,'k--','LineWidth',0.5);
    xline(zeroEv,'k--','LineWidth',0.5);
    title([fns{i}(1:end-5) '  %VE = ' num2str(round(varExp*100,2))])
    xlim([rez.time(1)+0.2,rez.time(end)])
    ax = gca;
    ax.FontSize = 20;
    ax.YTick = [];
    if i < 8
        ax.XTick = [];
    end
    ax.Color = [237, 237, 237]./255;
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





