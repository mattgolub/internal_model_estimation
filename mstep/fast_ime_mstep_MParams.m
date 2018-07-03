function estParams = fast_ime_mstep_MParams(E_X_posterior, COV_X_posterior, const)
% estParams = fast_ime_mstep_MParams(E_X_posterior, COV_X_posterior, const)
%
% @ Matt Golub, 2014.

[EXDim,T] = size(E_X_posterior);
xDim = const.xDim;
uDim = const.uDim;

%% Determine number of relevant timesteps in each trial for internal state estimates
num_internal_states_per_t = EXDim/xDim - 1;
N_x = num_internal_states_per_t*T;

%% Sums across state chains
sum_COV_X = sum(COV_X_posterior,3);

sumk_cov_xkxk = zeros(xDim);
sumk_cov_xkxkm1 = zeros(xDim);
sumk_cov_xkm1xkm1 = zeros(xDim);
for i = 1:xDim:(EXDim-xDim)
    idx_xkm1 = i:(i+xDim-1);
    idx_xk = idx_xkm1 + xDim;
    sumk_cov_xkxkm1 = sumk_cov_xkxkm1 + sum_COV_X(idx_xk,idx_xkm1);
    sumk_cov_xkm1xkm1 = sumk_cov_xkm1xkm1 + sum_COV_X(idx_xkm1,idx_xkm1);
    sumk_cov_xkxk = sumk_cov_xkxk + sum_COV_X(idx_xk, idx_xk);
end

sumk_cov_akak = [sumk_cov_xkm1xkm1 zeros(xDim,uDim+1);
    zeros(uDim,xDim+uDim+1);
    zeros(1,uDim+xDim+1)];

%% Gather sufficient statistics from E-step
E_xk = zeros(xDim,T*num_internal_states_per_t);
E_xkm1 = zeros(xDim,T*num_internal_states_per_t);
j = 1;
for i = 1:xDim:(EXDim-xDim)
    idx_xkm1 = i:(i+xDim-1);
    idx_xk = idx_xkm1 + xDim;
    t_idx = (T*(j-1)+1):(T*j);
    E_xk(:,t_idx) = E_X_posterior(idx_xk,:);
    E_xkm1(:,t_idx) = E_X_posterior(idx_xkm1,:);
    j = j + 1;
end
E_a = [E_xkm1; const.Urep; ones(1,T*num_internal_states_per_t)];
sum_Ex_Ea = E_xk*E_a';
temp = [E_xkm1*const.Urep' sum(E_xkm1,2)];
sum_Ea_Ea = [E_xkm1*E_xkm1' temp;
    temp' const.sum_U1_U1];

%% M variables update (AND V)
sum_E_aa = sumk_cov_akak + sum_Ea_Ea;
sum_E_xa = [sumk_cov_xkxkm1 zeros(xDim,uDim+1)] + sum_Ex_Ea;

% Update for M1, M2 and m0
% x_t = M*a_t, a_t = [x_{t-1}; u{t}; 1]
M = sum_E_xa / sum_E_aa;
estParams.M1 = M(:,1:xDim);
estParams.M2 = M(:,xDim+1:uDim+xDim);
estParams.m0 = M(:,end);

% Update for V
estParams.V = (1/N_x)*(diag(sumk_cov_xkxk) + diagProduct(E_xk,E_xk') - diagProduct(M,sum_E_xa'));