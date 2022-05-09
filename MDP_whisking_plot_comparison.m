function MDP_whisking_plot_comparison(o,labels)
rng default

q = [pi/4,pi/8,-pi/8,-pi/4]'; % Whisking angles

if ~iscell(o)
    o = {o};
end

for i = 1:numel(o)   
    y = q(o{i});
    Y = interp(y,8,2);
    Y = Y(1:end-8);
    X = [cos(Y),sin(Y)];
    X = X + spm_conv(randn(size(X)),2)/8;
    for t = 1:size(X,1)
        subplot(numel(o),2,i)
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

        subplot(2,2,2+i)
        plot(1:t,Y(1:t),'k')

        axis([1 size(X,1) -pi/2 pi/2]), axis square
        xlabel('Time (arbitrary units)')

        ylabel('Whisking angle')
        title('Whisking dynamics')

        pause(0.01)
        drawnow
    end

    subplot(numel(o),2,i)
    histogram2([X(:,1);-X(:,1)],[X(:,2);X(:,2)],-1.5:0.1:1.5,-1.5:0.1:1.5,'DisplayStyle','tile','ShowEmptyBins','on');

    axis square
    colormap gray

    title([labels{i}])
end

return