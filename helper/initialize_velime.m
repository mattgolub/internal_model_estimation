% See comments in velime_fit.m for descriptions of these settings.
TOL = 1e-8;
ZERO_TOL = 1e-10;
MAX_ITERS = 5e3;
MIN_ITERS = 10;
VERBOSE = false;
INIT_METHOD = 'current_regression';
INIT_PARAMS = NaN;
DO_LEARN_M_PARAMS = true;
DO_LEARN_ALPHA_PARAMS = true;
assignopts(who, varargin);

DT = 1; % DO NOT CHANGE!
% Without loss of generality, dt was assumed to be 1 throughout the 
% implementation. In some places (e.g., the e-step), the explicit value of 
% dt is used. In others (e.g., the m-step), it is omitted. If dt==1, 
% these inconsistent coding conventions still result in a self-consistent 
% model and fitting procedure. Otherwise, the model and fitting procedure 
% will no longer be self-consistent, and very bad things will happen.

% Initialization of parameters
% For {M1,M2,m0}
switch lower(INIT_METHOD)
    case 'target_regression'
        % Regress target direction against spikes.
        [B, b0, W_v] = target_regression(data);
        estParams.A = eye(2);
        estParams.B = B;
        estParams.b0 = b0;
        estParams.W_v = W_v;
        estParams.W_p = 0*W_v;
        estParams.TAU = TAU;
        estParams.dt = DT;
        % alphas initialized below
    case 'current_regression'
        [A, B, b0, W_v] = current_regression(data);
        % Regress current cursor-to-target direction against spikes.
        estParams.A = 0*A;
        estParams.B = B;
        estParams.b0 = b0;
        estParams.W_v = W_v;
        estParams.W_p = 0*W_v;
        estParams.TAU = TAU;
        estParams.dt = DT;
        % alphas initialized below
    case 'init_params'
        estParams = INIT_PARAMS;
        estParams.TAU = TAU;
        estParams.dt = DT;
    otherwise
        error('Unsupported initialization method');
end

[C,G,const] = velime_assemble_data(data,TAU,DT);

if ~isfield(estParams,'alpha')
    %% Initialize alpha's using prior distribution over latents
    [E_X, COV_X] = velime_prior_expectation(C, estParams);
    init_alphaParams = velime_mstep_alphaParams(G, E_X, repmat(COV_X,[1 1 const.T]), const);
    
    %% Assign estParams
    estParams.alpha = init_alphaParams.alpha;
    estParams.R = init_alphaParams.R;
end