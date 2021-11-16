function MDP = sensorimotor_exploration

rng('shuffle')

% Set Up:
T = 16;

N = 32;

% Specify model:
mdp   = generate_mdp_darkroom_1611(T);

% Specify number of trials:
M(1:N) = deal(mdp);
MDP = spm_MDP_VB_X(M);

s = MDP.s;
for t = 1:T
    if s(1,t) == 5
        s(1,t) = 3;
    end
    if s(1,t) == 6
        s(1,t) = 2;
    end
end
x = s(1,:)';
histfit(x)
title('Whisking positions');
xlim([0 5]);
xticklabels({'','Fully Protracted', 'Partially Protracted', 'Partially Retracted', 'Fully Retracted'});
ylim([0 10]);
ylabel('No. of Hidden State Visited')



spm_figure('GetWin','Figure1'); clf
spm_MDP_VB_trial(MDP(1));
