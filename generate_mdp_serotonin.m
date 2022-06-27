function mdp = generate_mdp_serotonin(T,zeta,rho,design,starting_position)

addpath 'D:\PhD\spm12'
addpath 'D:\PhD\spm12\toolbox\DEM'
addpath 'D:\PhD\model development\script'

%{

Description:

Free whisking is always under complete uncertainty (dark room problem).
When next to the platform, fully and partially protracted positions generate
2 different sensory outcomes, which are used for identifying the 
object. These two protracted positions perceive the sensory
outcomes under different probability distribution, where the partially protracted hidden
state is fully precise.
%}
% Set up:

% beta = 0.10;
s(2,1) = design; % Experimental design whether the agent is next to (1) or far away (2)

% Initial matrix D P(s)

D{1} = starting_position; %[protraction2; protraction1;
% retraction 1; retraction2; retraction1; protraction1]

D{2} = [0 1]'; % Context[next to, or far from, the object]

% Hidden factors and states:
Nf = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
% Outcomes and outcome modalities:
No    = [3 4];
% Description of No respectively:
%{
No(1) = [Nothing sensed; Edge; Surface]
No(2) = [fully protracted, partially protracted, partially retracted, partially protracted]
%}
Ng    = numel(No);

% Likelihood matrix A P(o_t|s_t)

% Initialising likelihood:
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
    a{g} = A{g};
end
for f = 1:Nf
    B{f} = zeros(Ns(f));
end

for f1 = 1:Ns(1) % Body control [completely protracted <-> completely retracted]
    for f2 = 1:Ns(2) % Context [Next to or far away from the object]

        % Generative process
        A{1}(:,f1,2) = 1/3;
        A{1}(2:3,1,1) = 0.5;
        A{1}(2,[2 6],1) = 1;
        A{1}(:,3:5,1) = 1/3;
        % Protraction/Retraction
        A{2}(:,1:4,f2) = eye(No(2));
        A{2}(2,6,f2) = 1;
        A{2}(3,5,f2) = 1;

        % Generative model
        a{1}(:,f1,2) = spm_softmax(0*log(A{1}(:,f1,2)+exp(-8)));
        % partially protracted hidden state
        a{1}(:,:,1) = spm_softmax(zeta*log(A{1}(:,:,1)+exp(-8)));
        % fully protracted hidden state
        %         a{1}(:,[1 4],1) = spm_softmax(0.5*log(A{1}(:,[1 4],1)+exp(-8)));

        a{2}(:,f1,f2) = A{2}(:,f1,f2);

        % Transition matrix P(s_t+1|s_t,pi)
        B{1}(:,:,1) = circshift(eye(Ns(1)),1);
        B{1}(2,[1 3 5],2) = 1;
        B{1}(5,[2 4 6],2) = 1;

        B{2} = eye(Ns(2));
    end
end
% b{1} = B{1};
% b{2} = spm_softmax(rho*log(B{2} +exp(-8)));

E = [0.6 0.4]';
E = spm_softmax(rho*log(E));
%% Specify 1-step policies
U(:,:,1) = [1 2];
U(:,:,2) = [1 1];


%% MDP structure
%=================================================
mdp.T = T;
mdp.U = U;
mdp.A = A;
mdp.a = a;
mdp.B = B;
% mdp.b = b;
mdp.D = D;
mdp.E = E;
% mdp.e = e;
mdp.s = s;
mdp.eta = 0;
mdp.beta = 1;


%% Labels
label.factor{1} = 'Whisking Cycles';
label.name{1} = {'Fully Protracted','Partially Protracted','Partially Retracted','Fully Retracted','Partially Retracted','Partially Protracted'};
label.factor{2} = 'Context';
label.name{2} = {'Near','Far'};

label.modality{1} = 'Sensory Information';
label.outcome{1} = {'Nothing Sensed','Edge','Surface'};
label.modality{2} = 'Protraction/Retraction';
label.outcome{2} = {'Fully Protracted','Partially Protracted','Partially Retracted','Fully Retracted'};
label.action{1} = {'Large','Small'};
mdp.label = label;

return