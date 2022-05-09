function MDP = whiskers_main

rng('default')

% Set Up:
N = 1; % All simulations apart timeplot shows only the first simulation trial
T = 16;
design = 1; % 1 == near to the platform; 2 == far from the platform

Sim_habits     = 0; % Simulate a modulation of habits if 1
% Adjust ZETA to have only per one entry; RHO has 3 entries
Sim_precision  = 0; % Simulate a modulation of somatosensory precision if 1
% Adjust RHO to have only per one entry; ZETA has 3 entries
Sim_comparison = 0; % Simulate a face validity of far and near conditions if 1
% Adjust ZETA and RHO to have only one entry
Sim_timeplot   = 1; % Simulate a range of all precision interactions if 1

% Always simulate only one condition and adjust precision terms adequately

ZETA  = [0.1 0.5 2];
nze   = numel(ZETA);
RHO = [0.1 0.5 2];
nrho = numel(RHO);
store_time = zeros(nze,nrho);

% Face validity

if Sim_comparison == 1

   design = 1;
   zeta = 0.5;
   rho = 0.5;
    mdp = generate_mdp_serotonin(T,zeta,rho,design);
    MDP = spm_MDP_VB_X(mdp);
    o{1} = MDP.o(2,:);
    
    design = 2;
    mdp = generate_mdp_serotonin(T,zeta,rho,design);
    MDP = spm_MDP_VB_X(mdp);
    o{2} = MDP.o(2,:);

    labels = {'Near','Far'};
    spm_figure('GetWin','Whisking plot'); clf

    MDP_whisking_plot_comparison(o,labels)
end
for ze = 1:nze
    for e = 1:nrho
        zeta  = ZETA(ze);
        rho = RHO(e);

        mdp = generate_mdp_serotonin(T,zeta,rho,design);

        M(1:N) = deal(mdp);
        MDP = spm_MDP_VB_X(M);
        % Hypothesis 1: Serotonin modulates tactile precision
        if Sim_precision == 1 && nrho == 1
            o{ze} = MDP(1).o(2,:);
            if ze == 1
                spm_figure('GetWin','LFPs for tactile precision = 0.1')
                spm_MDP_VB_LFP(MDP(1));
            elseif ze == 2

                spm_figure('GetWin','LFPs for tactile precision = 0.5')
                spm_MDP_VB_LFP(MDP(1));
            elseif ze == 3
                spm_figure('GetWin','LFPs for tactile precision = 2')
                spm_MDP_VB_LFP(MDP(1));
            end
            labels = {'Low TP','Medium TP','High TP'};

        end
        % Hypothesis 2: Serotonin modulates the prior policy precision
        if Sim_habits == 1 && nze == 1
            o{e} = MDP.o(2,:);
            if e == 1
                spm_figure('GetWin','LFPs for Habitual prior precision = 0.1')
                spm_MDP_VB_LFP(MDP);
            elseif e == 2

                spm_figure('GetWin','LFPs for Habitual prior precision = 0.5')
                spm_MDP_VB_LFP(MDP);
            elseif e == 3
                spm_figure('GetWin','LFPs for Habitual prior precision = 2')
                spm_MDP_VB_LFP(MDP);
            end
            labels = {'Low Habits','Medium Habits','Strong Habits'};

        end

        if Sim_timeplot == 1 && design == 1
            for n = 1:N
                t = 1;
                while t < T - 2 && not(MDP(n).u(1,t) == 2 && MDP(n).u(1,t+2) == 2)
                    t = t + 1;
                end

                if t < T - 2 && MDP(n).u(1,t) == 2 && MDP(n).u(1,t+2) == 2
                    store_time(ze,e) = store_time(ze,e) + t;
                else
                    t = 16;
                    store_time(ze,e) = store_time(ze,e) + t;
                end

            end
            store_time(ze,e) = store_time(ze,e)./N;

        elseif Sim_timeplot == 1 && design == 2
          
            for n = 1:N
                t = 2;
                while t < T && MDP(n).u(1,t-1) == 1 && MDP(n).u(1,t) == 1
                    t = t + 1;
                end

                if t < T && MDP(n).u(1,t) == 2
                    store_time(ze,e) = store_time(ze,e) + t;
                else
                    t = 16;
                    store_time(ze,e) = store_time(ze,e) + t;
                end

            end
            store_time(ze,e) = store_time(ze,e)./N;
        end

    end
end

if Sim_timeplot == 1

    bar(store_time(:,:), 'grouped');


    title('Time of a whisking switch');
    xlim([0 4]);
    xticklabels({ZETA});
    xlabel('Tactile precision');
    ylim([0 19]);
    ylabel('Average time')
    legend('HPP = 0.1', 'HPP = 0.5', 'HPP = 2');

elseif Sim_timeplot == 0 && Sim_comparison == 0

    spm_figure('GetWin','Whisking plot'); clf

    MDP_whisking_plot(o,labels)
end