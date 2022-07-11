function [u,v] = LFP_serotonin_tactile(MDP,ze,T)

% This function is taken from spm_MDP_VB_LFP and the original copright
% belongs to Karl Friston. See the original for more details.

% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_LFP.m 7653 2019-08-09 09:56:25Z karl $

f = 1;
UNITS = [];
SPECTRAL = 0;


% dimensions
%--------------------------------------------------------------------------
Nt     = length(MDP);               % number of trials
try
    Ne = size(MDP(1).xn{f},4);      % number of epochs
    Nx = size(MDP(1).B{f}, 1);      % number of states
    Nb = size(MDP(1).xn{f},1);      % number of time bins per epochs
catch
    Ne = size(MDP(1).xn,4);         % number of epochs
    Nx = size(MDP(1).A, 2);         % number of states
    Nb = size(MDP(1).xn,1);         % number of time bins per epochs
end

% units to plot
%--------------------------------------------------------------------------
ALL   = [];
for i = 1:Ne
    for j = 1:Nx
        ALL(:,end + 1) = [j;i];
    end
end
if isempty(UNITS)
    UNITS = ALL;
end

% summary statistics
%==========================================================================
for i = 1:Nt

    % all units
    %----------------------------------------------------------------------
    str    = {};
    try
        xn = MDP(i).xn{f};
    catch
        xn = MDP(i).xn;
    end
    for j = 1:size(ALL,2)
        for k = 1:Ne
            zj{k,j} = xn(:,ALL(1,j),ALL(2,j),k);
            xj{k,j} = gradient(zj{k,j}')';
        end
        str{j} = sprintf('%s: t=%i',MDP(1).label.name{f}{ALL(1,j)},ALL(2,j));
    end
    z{i,1} = zj;
    x{i,1} = xj;

    % selected units
    %----------------------------------------------------------------------
    for j = 1:size(UNITS,2)
        for k = 1:Ne
            vj{k,j} = xn(:,UNITS(1,j),UNITS(2,j),k);
            uj{k,j} = gradient(vj{k,j}')';
        end
    end
    v{i,1} = vj;
    u{i,1} = uj;
    %             str{j} = sprintf('%s: t=%i',MDP(1).label.name{f}{ALL(1,j)},ALL(2,j));
end
z{i,1} = zj;
x{i,1} = xj;
dt  = 1/64;                              % time bin (seconds)
t   = (1:(Nb*Ne*Nt))*dt;                 % time (seconds)
Hz  = 4:32;                              % frequency range
n   = 1/(4*dt);                          % window length
w   = Hz*(dt*n);
LFP = spm_cat(x);
wft = spm_wft(LFP,w,n);
csd = sum(abs(wft),3);
lfp = sum(LFP,2);
phi = spm_iwft(sum(wft(1,:,:),3),w(1),n);
lfp = 4*lfp/std(lfp) + 16;
phi = 4*phi/std(phi) + 16;

hold on
subplot(3,3,6+ze)
imagesc(t,Hz,csd), axis xy
hold on
plot(t,phi,'w')
hold off

ax = gca;
ax.XTickLabel = ax.XTick*1000/T*1.2;

ax.YTickLabel = ax.YTick*0.1;
caxis([0 0.1])

lgnd = legend('\color{white} LFP','box', 'off',Location='northeast');
set(lgnd,'color','none')
lgnd.ItemTokenSize = [10,10];
xlabel('Time (ms)','fontsize',12), ylabel('frequency (Hz)','fontsize',12)
if Nt == 1, axis square, end