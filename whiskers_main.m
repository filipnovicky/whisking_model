function MDP = whiskers_main

clear; clc

% Add a path to the spm12 file (see an example):
addpath 'D:\PhD\spm12'
addpath 'D:\PhD\spm12\toolbox\DEM'
addpath 'D:\PhD\model development\script'

% Set a file where you want to save the data (see an example):
path =  'D:\PhD\model development\script\results';

if ~exist(path, 'dir')
    mkdir(path)
end

rng('default')
% Select a simulation
% RUN ONLY ONE SIMULATION PER TIME

Sim_habits     = 0; % Simulate a modulation of habits if 1
Sim_precision  = 0; % Simulate a modulation of somatosensory precision if 1
Sim_comparison = 0; % Simulate a face validity of far and near conditions if 1
Sim_timeplot   = 0; % Simulate a range of all precision interactions if 1


if Sim_habits == 1
    % Do not change
    Sim_timeplot = 0;
    Sim_precision = 0;
    Sim_comparison = 0;
    N = 1;
    T = 16;
    % Feel free to change:
    design = 1; % 1 == near condition; 2 == far condition
    RHO = [0.01 1 2];
    ZETA = 0.5;

elseif Sim_precision == 1
    % Do not change:

    Sim_timeplot = 0;
    Sim_habits = 0;
    Sim_comparison = 0;
    N = 1;
    T = 16;
    % Feel free to change:
    design = 2; % 1 == near condition; 2 == far condition
    RHO = 0.5;
    ZETA = [0.01 0.5 2];

elseif Sim_timeplot == 1
    % Do not change:
    Sim_precision = 0;
    Sim_habits = 0;
    Sim_comparison = 0;
    % Feel free to change:
    design = 2; % 1 == near condition; 2 == far condition
    ZETA = fliplr([0.01 0.1 0.3 0.5 1 1.5 2]);
    RHO = [0.01 0.1 0.3 0.5 1 1.5 2 5 10 100];
    N = 1;
    T = 16;
end

% Random whisker protraction position for temporal analysis
starting_position = [1 0 0 0 0 0]';
if Sim_timeplot == 1
    starting_position = [1/3 1/3 0 0 0 1/3]';
end


% Face validity

if Sim_comparison == 1

    %Set up:
    T = 16;
    N = 1;
    zeta = 1;
    rho = 1.5;
    %Set-up for the far condition:
    design = 1;


    mdp = generate_mdp_serotonin(T,zeta,rho,design,starting_position);
    MDP = spm_MDP_VB_X(mdp);
    o{1} = MDP.o(2,:);
    spm_figure('GetWin','Whisking plot'); clf
    subplot(3,2,5)
    policy_posterior = MDP.R;
    xq = 0:0.01:15;

    for i = 1:2
        hold on
        plot(interp1(policy_posterior(i,:),xq,"spline"))
    end
    box on
    legend('Small Amplitude', 'Large Amplitude')

    title('Posterior Policy Probability', 'fontsize',12)

    xlabel('Time (ms)','fontsize',12)
    ax = gca;
    ylim([0 1.5])
    ax.XTickLabel = ax.XTick*0.3125*T*0.04; % Set up the time on the x-axis

    % Set-up for the far condition:
    design = 2;


    mdp = generate_mdp_serotonin(T,zeta,rho,design,starting_position);
    MDP = spm_MDP_VB_X(mdp);
    o{2} = MDP.o(2,:);
    hold on
    subplot(3,2,6)
    policy_posterior = MDP.R;
    for i = 1:2
        hold on
        plot(interp1(policy_posterior(i,:),xq,"spline"))
    end
    box on

    title('Posterior Policy Probability','fontsize',12)

    xlabel('Time (ms)','fontsize',12)
    legend('Small Amplitude', 'Large Amplitude')
    ylim([0 1.5])
    ax = gca;
    ax.XTickLabel = ax.XTick*0.3125*T*0.04; % Set up the time on the x-axis


    labels = {'Whisker Position','Whisker Position'};
    hold on

    MDP_whisking_plot_comparison(o,labels,T)
end

if Sim_comparison == 0

    nze   = numel(ZETA);
    nrho  = numel(RHO);
    store_time = zeros(nze,nrho,N);
    av_time = zeros(nze,nrho);
    median_time = zeros(nze,nrho);

    for ze = 1:nze
        for e = 1:nrho
            zeta  = ZETA(ze);
            rho = RHO(e);

            mdp = generate_mdp_serotonin(T,zeta,rho,design,starting_position);

            M(1:N) = deal(mdp);
            MDP = spm_MDP_VB_X(M);
            % Hypothesis 1: Serotonin modulates tactile precision
            if Sim_precision == 1 && nrho == 1
                o{ze} = MDP(1).o(2,:);
                % Simulate LFP signals for pre-specified zetas
                if ze == 1
                    spm_figure('GetWin','LFP figure')
                    LFP_serotonin_tactile(MDP(1),ze,T)
                elseif ze == 2
                    hold on
                    LFP_serotonin_tactile(MDP(1),ze,T)

                elseif ze == 3
                    hold on
                    LFP_serotonin_tactile(MDP(1),ze,T)

                end
                labels = {'Low tactile precision','Medium tactile precision','High tactile precision'};

            end
            % Hypothesis 2: Serotonin modulates the prior policy precision
            if Sim_habits == 1 && nze == 1
                o{e} = MDP.o(2,:);
                % Simulate LFP signals for pre-specified rhos
                if e == 1
                    spm_figure('GetWin','LFP figure')

                    LFP_serotonin_tactile(MDP(1),e,T)

                elseif e == 2
                    hold on
                    LFP_serotonin_tactile(MDP(1),e,T)

                elseif e == 3
                    hold on
                    LFP_serotonin_tactile(MDP(1),e,T)

                end
                labels = {'Weak habitual precision','Medium habitual precision','Strong habitual precision'};
            end

            % Store data to simulate heatmaps
            small_whisk = false; 

            if Sim_timeplot == 1 && small_whisk == false
                for n = 1:N
                    t = 1;
                    while t < T - 2 && not(MDP(n).u(1,t) == 2 && MDP(n).u(1,t+2) == 2)
                        t = t + 1;
                    end

                    if t < T - 2 && MDP(n).u(1,t) == 2 && MDP(n).u(1,t+2) == 2
                        store_time(ze,e,n) = store_time(ze,e,n) + t;
% Save the time of the switch when two small amplitude actions are selected
                    else
                        t = 16; % save 16 instead if no switch occurred
                        store_time(ze,e,n) = store_time(ze,e,n) + t;
                    end

                end

            elseif Sim_timeplot == 1 && small_whisk == true && design == 2
 % This considers a switch only after one small whisk cycle in the far
 % condition
                
                                for n = 1:N
                                    t = 2;
                                    while t < T && MDP(n).u(1,t-1) == 1 && MDP(n).u(1,t) == 1
                                        t = t + 1;
                                    end
                
                                    if t < T && MDP(n).u(1,t) == 2
                                        store_time(ze,e,n) = store_time(ze,e,n) + t;
 % Save the time of the switch when two small amplitude actions are selected
                                    else
                                        t = 16;
                                        % save 16 instead if no switch occurred
                                        store_time(ze,e,n) = store_time(ze,e,n) + t;
                                    end
                
                                end
            end
            av_time(ze,e) = mean(store_time(ze,e,:));
            median_time(ze,e) = median(store_time(ze,e,:));
        end

    end

end

if Sim_timeplot == 1 % Simulate discrete heatmaps

    yvalues = ZETA;
    xvalues = RHO;
    h = heatmap(xvalues,yvalues,median_time,'CellLabelColor','none','Colormap',parula,'ColorLimits',[0 16],'fontsize',12);
    h.Title = 'Median time delay of the behavioural switch';
    h.YLabel = 'Tactile precision';
    h.XLabel = 'Habitual precision';
    annotation('textarrow',[0.95,0.95],[0.5,0.5],'string','Time Delay','fontsize',12, ...
        'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90);
    hHeatmap = struct(h).Heatmap;
    hHeatmap.GridLineStyle = 'none';

    figure
    h = heatmap(xvalues,yvalues,av_time,'CellLabelColor','none','Colormap',parula,'ColorLimits',[0 16],'fontsize',12);
    h.Title = 'Average time delay of the behavioural switch';
    h.YLabel = 'Tactile precision';
    h.XLabel = 'Habitual precision';
    annotation('textarrow',[0.95,0.95],[0.5,0.5],'string','Time Delay','fontsize',12, ...
        'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90);
    hHeatmap = struct(h).Heatmap;
    hHeatmap.GridLineStyle = 'none';

elseif Sim_timeplot == 0 && Sim_comparison == 0
    hold on
    MDP_whisking_plot(o,labels,T)

end

% Save the data for heatmaps (see an example):
if Sim_timeplot == 1
if design == 1
    save(strcat(path,'\near\av_time.mat'),'av_time');
    save(strcat(path,'\near\median_time.mat'),'median_time');
elseif design == 2
    save(strcat(path,'\far\av_time.mat'),'av_time');
    save(strcat(path,'\far\median_time.mat'),'median_time');
end
end