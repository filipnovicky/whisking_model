function MDP_whisking_plot(o,labels,T)

% This function takes a vector of whisking outcomes (i.e., proprioceptive
% outcomes associating with whisking position). These are indexed
% numerically such that o = 1 implies full protraction and o = 4 implies
% full retraction. If o is a cell array, each element is plotted in an
% alternative subplot.
%--------------------------------------------------------------------------

rng default% for reproducibility

q = [pi/4,pi/8,-pi/8,-pi/4]'; % Whisking angles

if ~iscell(o)
    o = {o};
end

for i = 1:numel(o)
    y = q(o{i});
    Y = interp(y, 8, 2);
    Y = Y(1:end-8);
    X = [cos(Y), sin(Y)];
    X = X + spm_conv(randn(size(X)),2)/8;
    for t = 1:size(X,1)
        subplot(numel(o),3,i)
        if t > 3
            plot([0 X(t-2,1)],[0 X(t-2,2)],'Color',[0.9 0.9 0.9]), hold on
            plot([0 -X(t-2,1)],[0 X(t-2,2)],'Color',[0.9 0.9 0.9]), hold on
        end
        if t > 3
            plot([0 X(t-1,1)],[0 X(t-1,2)],'Color',[0.5 0.5 0.5]), hold on
            plot([0 -X(t-1,1)],[0 X(t-1,2)],'Color',[0.5 0.5 0.5]), hold on

        end

        plot([0 X(t,1)],[0 X(t,2)],'Color',[0.1 0.1 0.1]), hold on
        plot([0 -X(t,1)],[0 X(t,2)],'Color',[0.1 0.1 0.1]), hold off
       
        axis([-1.5 1.5 -1.5 1.5]), axis square
        title([labels{i}])
        subplot(3, 3, 3+i)
        plot(1:t,Y(1:t),'k')
        axis([1 size(X,1) -pi/2 pi/2]), axis square
        xlabel('Time (ms)','fontsize',12)
        ax = gca;
        ax.XTickLabel = ax.XTick * 0.3125*T/2;
        ylabel('Whisking angle','fontsize',12)
        pause(0.01)
        drawnow

    end
    subplot(numel(o),3,i)
    histogram2([X(:,1);-X(:,1)],[X(:,2);X(:,2)],-1.5:0.1:1.5,-1.5:0.1:1.5,'DisplayStyle','tile','ShowEmptyBins','on');
    axis square
    ylabel('Whisking angle','fontsize',12)
    colormap gray
    title([labels{i}])

end


return
