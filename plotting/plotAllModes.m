function plotAllModes(rez,ev,alignEv,plt)

% get field names for each mode
[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');

psth = rez.psth - mean(rez.psth,1);

sample = mode(ev.sample) - mode(ev.(alignEv));
delay  = mode(ev.delay) - mode(ev.(alignEv));
zeroEv  = 0; % corresponds to either gocue or moveonset

% plot each mode
fig = figure;
for i = 1:numel(fns)
    subplot(5,2,i);
    hold on
    for j = 1:numel(plt.conditions)
        cond = plt.conditions(j);
        latent = mySmooth(psth(:,:,cond)*rez.(fns{i}),plt.smooth);
        plot(rez.time,latent,'Color',plt.colors{j},'LineWidth',plt.lw(j));
    end
%     xline(sample,'k--','LineWidth',0.5);
%     xline(delay,'k--','LineWidth',0.5);
    xline(zeroEv,'k--','LineWidth',0.5);
    title(fns{i}(1:end-5),'FontSize',13)
    xlim([-2.5,2.1])
    ax=gca;
    ax.FontSize=11;
    hold off
end

leg = legend(plt.legend);
leg.Position = [0.60,0.066,0.28,0.13];
leg.FontSize = 14;
sgtitle(plt.title)

h=axes(fig,'visible','off');
h.XLabel.Visible='on';
h.YLabel.Visible='on';
ylabel(h,'Activity (a.u.)','FontSize',14);

xlabel(h,['Time (s) from ' alignEv],'FontSize',14);

end % plotAllModes





