function [C,G,const] = velime_assemble_data(data, TAU, dt)
% [C,G,const] = velime_assemble_data(U, P, TARGETS, TAU, dt)
%
% Produce structured data matrices required for the EM algorithm.
%
% INPUTS: 
%
%
% OUTPUTS:
% C: has columns c_t = [p_{t-tau}; v_{t-tau-1}; u_{t-tau+1}; ... ; u_{t+1}]
%       ( where v_{t-tau-1} = (p_{t-tau} - p_{t-tau-1})/dt )
%
% G: has columns g_t = target position from appropriate trial

U = data.spike_counts;
P = data.cursor_position;
TARGETS = data.target_position;

t_start = TAU+2; % "whiskers" are only defined for x_k^t with t>=TAU+2
num_trials = numel(U);
T = cellfun(@(u)(size(u,2)),U);

xDim = 4;
gDim = 2; 
uDim = mode(cellfun(@(u)(size(u,1)),U));
EXDim = (TAU+1)*xDim;

% Write this, give xMap: one string for each set of indices into x
[K, x_idx] = velime_x_index(TAU,xDim);

T_valid = sum(max(0,T-(TAU+2)));
C = zeros(xDim + (TAU+1)*uDim,T_valid);
G = zeros(gDim,T_valid);

idx = 1;
for trialNo = 1:num_trials
    const.trial_map{trialNo} = idx:((idx-1)+T(trialNo)-t_start);
    V_trial = diff(P{trialNo},1,2)/dt; % velocities computed from differences in positions
    for t = t_start:T(trialNo)-1
        U_seq = U{trialNo}(:,t-TAU+1:t+1);
        
        % Available from feedback: p(t-tau) and v(t-tau-1) = p(t-tau)-p(t-tau-1)
        %   (using convention: p(t) = p(t-1) + v(t-1)dt)
        % Sanity check, these should be identical:
        %   V_trial(:,t-TAU-1)
        %   (P{trialNo}(:,t-TAU)-P{trialNo}(:,t-TAU-1))/dt
        
        C(:,idx) = [P{trialNo}(:,t-TAU); V_trial(:,t-TAU-1); U_seq(:)];
        G(:,idx) = TARGETS{trialNo}(:,1);
        idx = idx + 1;
    end
end

const.T = T_valid;
const.dt = dt;
const.uDim = uDim;
const.gDim = gDim;
const.xDim = xDim;
const.pDim = 2;
const.vDim = 2;
const.EXDim = EXDim;
const.K = K;
const.x_idx = x_idx;
const.x_pt_idx = EXDim-3:EXDim-2;
const.x_vt_idx = EXDim-1:EXDim;

% Here TAU+1 gives the number of timesteps worth of u's that contribute to
% each latent state chain, x.
const.Urep = zeros(uDim,T_valid*(TAU+1));
for j = 1:TAU+1
    % j=1, xDim=4 --> 5:(4+uDim)
    % j=2, xDim=4 --> (5+uDim):(4+2*uDim)
    idx_u = (xDim + 1 + uDim*(j-1)):(xDim + uDim*j); % rows
    idx_t = (T_valid*(j-1)+1):T_valid*j; % columns
    const.Urep(:,idx_t) = C(idx_u,:);
end

sum_Urep = sum(const.Urep,2);
const.sum_U1_U1 =  [const.Urep*const.Urep' sum_Urep;
    sum_Urep' T_valid*(TAU+1)];