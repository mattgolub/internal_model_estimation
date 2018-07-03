function sweep = velime_sweep_TAU(data, TAUs, varargin)
% sweep = velime_sweep_TAU(data, TAUs, ...)
%
% INPUTS:
%
%   data.spike_counts: 
%       Binned spike counts. [1 x # trials] cell array. Each element is
%       [# neurons x # timesteps].
%
%   data.cursor_position: 
%       2D cursor positions. [1 x # trials] cell array. Each element is
%       [2 x # timesteps].
%
%   data.target_position: 
%       2D target positions. [1 x # trials] cell array. Each element is 
%       [2x1] or [2 x # timesteps].
%
%   TAUs: 
%       A vector containing the set of candidate TAU settings. IME models 
%       are fit and evaluated for each value of TAU in TAUs. TAU is 
%       the sensory feedback delay, measured in timesteps. Each element
%       must be a non-negative integer.
%
% OUTPUTS: 
%   sweep: 
%       A struct array, where element i corresponds to an IME model (or 
%       multiple models if using cross validation) with TAU = TAUs(i). Each
%       element contains the following fields:
%           estParams: 
%               Struct containing fit IME parameters. If using cross
%               validation (default is no), estParams is a struct array
%               where each element contains IME parameters fit to a 
%               different fold of training data. See velime_fit for 
%               descriptions of the IME parameters.
%           angular_error:
%               Struct with fields .model and .cursor, each of which are
%               [1 x # trials] vectors containing angular errors 
%               according to the extracted internal models. See
%               velime_evaluate for additional details.
%           LL:
%               Scalar value indicating the log-likelihood of the data
%               according to the extracted internal model(s). See optional
%               input DO_COMPUTE_LL below and velime_evaluate for more
%               details.
%           TAU:
%               The value of TAU used when fitting and evaluating the
%               corresponding IME models. This value is simply copied from
%               the corresponding entry in TAUs
%
% OPTIONAL INPUTS (and default values):
%
%   DO_CROSS_VALIDATE (false):
%       Logical indicating whether or not to perform cross validation
%       (i.e., report angular errors and LL based on trials held out during
%       model fitting (see velime_cross_validate for details of the cross
%       validation procedure). For the purposes of model selection across
%       settings of TAU, cross validation is not needed to control for
%       overfitting because the number of model parameters is not
%       determined by TAU.
%
%   DO_COMPUTE_LL (true):
%       Logical indicating whether or not to compute the data
%       log-likelihood. If setting DO_CROSS_VALIDATE to true, substantial 
%       additional computation is introduced to this procedure.
%
%   Additional optional arguments specific to velime_fit that are passed in
%   will be passed along to velime_fit (e.g., VERBOSE, MAX_ITERS).
%
% @ Matt Golub, 2018.

DO_CROSS_VALIDATE = false;
DO_COMPUTE_LL = true;
velime_args = assignopts(who,varargin);

T_START = 2+max(TAUs);

for TAU_idx = 1:numel(TAUs)
    TAU = TAUs(TAU_idx);
    
    fprintf('TAU = %d\n',TAU);
    
    if DO_CROSS_VALIDATE
        [estParams, ~, evaluations, ~] = velime_cross_validate(data, TAU, ...
            'T_START', T_START, ...
            'DO_COMPUTE_LL', true, ...
            velime_args{:});
        
        angular_error = evaluations.angular_error;
        LL = evaluations.LL;
    else
        estParams = velime_fit(data, TAU, velime_args{:});
        [E_P, E_V] = velime_predict(data, estParams);
        [angular_error, LL] = velime_evaluate(data, E_P, E_V, estParams, 'T_START', T_START, 'DO_COMPUTE_LL', DO_COMPUTE_LL);
    end
    
    sweep(TAU_idx).estParams = estParams;
    sweep(TAU_idx).angular_error = angular_error;
    sweep(TAU_idx).LL = LL;
    sweep(TAU_idx).TAU = TAU;
end