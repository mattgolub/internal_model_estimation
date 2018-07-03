function [angular_error, LL] = velime_evaluate(data, E_P, E_V, estParams, varargin)
% errors = velime_evaluate(data, E_P, E_V, estParams, ...)
% [errors, LL] = velime_evaluate(data, E_P, E_V, estParams, ...)
%
% Compute timestep-by-timestep angular errors in the actual cursor
% trajectory and according to the subject's internal model.
%
% INPUTS:
%
% data: struct with the following fields:
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
% E_P and E_V:
%   [1 x # trials] cell arrays containing the subject's internal estimates
%   of cursor position and velocity, respectively, as produced by
%   velime_extract_prior_whiskers. Each element is
%   [2*(TAU+1) x # timesteps].
%
% OUTPUTS:
%
% angular_error: struct with the following fields
%
%       angular_error.cursor: [1 x # trials] vector containing the average
%       absolute angular error within each trial. Here, angular error at
%       timestep t is defined as the angle by which the cursor would have
%       missed the target had it continued from position p(t) in the
%       direction of velocity v(t). The angular error is defined as zero if
%       velocity v(t) would have driven the cursor to visually overlap the
%       target. This is the black error annotation in the illustration of
%       Figure 3B in the paper.
%
%       angular_error.model: [1 x # trials] vector containing the average 
%       absolute errors within each trial, according to the internal model. 
%       Here, angular error at timestep t is defined as the angle by which 
%       the cursor would have missed the target had it continued from the
%       subject's internal estimate of cursor position \tilde{p}_t^t along
%       the direction of the subject's intended cursor velocity
%       \tilde{v}_t^t, both of which are based on the subject's internal
%       model. Zero errors are defined as in errors.cursor. This is the red
%       error annotation in the illustration of Figure 3B in the paper.
%
% LL (optional):
%       Scalar value indicating the log-likelihood of the data.
%
% OPTIONAL INPUTS (and default values):
%
% VERBOSE (false):
%       Logical indicating whether or not to print status updates to the
%       screen.
%
% T_START (TAU + 2):
%       Integer indicating the timestep to begin evaluating the data
%       (relevant for both angular errors and log-likelihood; the same
%       timestep is indicated for all trials in data). Must be >= TAU + 2
%       because "whiskers" are undefined for t<(TAU+2) (when some of the
%       requisite position and neural data do not exist, i.e., those data
%       are from timesteps before the trial began). Further, the subject
%       likely has not yet visually processed the target location, and thus
%       the assumption the subject is intending to drive the cursor
%       straight to the target is not valid.
%
%       When performing a sweep of TAU to determine the value that best
%       describes the data, set T_START to 2 + MAX_TAU, where MAX_TAU is
%       the maximum value of TAU being considered. This setting ensures
%       that all models (i.e., for each candidate TAU) are evaluated on
%       exactly the same data (i.e., exactly the same timesteps in each
%       trial). Note that the larger MAX_TAU is, the less data are
%       available, and the data that are available are more concentrated
%       toward the end of the trials (and thus may have different
%       characteristics than the data available with smaller values of
%       T_START).
%
% DO_COMPUTE_LL (false):
%       Logical indicating whether or not to compute the log-likelihood of
%       the data.
%
% DATA_ARE_TRAINING_DATA (false):
%       Logical indicating whether or not "estParams" were fit to "data."
%       This determines how the log-likelihood must be computed. When
%       properly cross-validating, estParams are fit to training data, and
%       evaluations are performed on held-out "test" or "validation" data.
%
% MAX_ITERS (50):
%       Integer indicating the number of EM iterations to perform if
%       necessary to obtain test-data-specific alpha(t) values. This is
%       only relevant for computing the log-likelihood of the data if
%       "estParams" was not fit to "data" (as indicated by 
%       DATA_ARE_TRAINING_DATA). In this case, EM iterations are performed
%       with a fixed set of internal model parameters, thus only learning
%       the auxiliary scale variables, alpha(t), and a corresponding
%       noise variance, R). Because the internal model parameters are not
%       learned, MAX_ITERS can be much smaller than is typically necessary
%       when learning all IME parameters.
%
% See also velime_fit, velime_predict, example1, example2, example3.
%
% @ Matt Golub, 2018

VERBOSE = false;
T_START = estParams.TAU + 2;
DO_COMPUTE_LL = false;
DATA_ARE_TRAINING_DATA = false;
MAX_ITERS = 50;
assignopts(who,varargin);

if estParams.TAU<0
    error('TAU must be >= 0');
end
if T_START < (estParams.TAU+2)
    error('T_START must be >= TAU + 2');
end

angular_error = compute_angular_errors(data, E_P, E_V, estParams, T_START);

if nargout == 2
    if DO_COMPUTE_LL
        % Compute the log-likelihood of the data
        LL = compute_LL(data, estParams, T_START, DATA_ARE_TRAINING_DATA, MAX_ITERS, VERBOSE);
    else
        LL = nan;
    end
end

end

function errors = compute_angular_errors(data, E_P, E_V, estParams, T_START)

N_DIMS = size(estParams.A,1);
TAU = estParams.TAU;
DT = estParams.dt;

acceptance_zone_radius = data.cursor_radius + data.target_radius;
% The "acceptance zone" is the circular region about the target center that 
% the cursor center must enter in a successful trial. The codepack's
% example data are from an experiment with cursor and target radii of 7mm 
% each. Trials were deemed successful if the cursor visibly overlapped the
% target (for 50ms). Thus acceptance_zone_radius = 14mm.

p_t_t_idx = ((TAU+1)*N_DIMS) + [-1 0];
v_t_t_idx = p_t_t_idx;
num_trials = numel(data.cursor_position);

errors.model = nan(1,num_trials);
errors.cursor = nan(1,num_trials);

for trial_idx = 1:num_trials
    G = data.target_position{trial_idx};
    
    T = size(data.cursor_position{trial_idx},2);
    P_t = data.cursor_position{trial_idx};
    P_tp1 = [data.cursor_position{trial_idx}(:,2:end) nan(N_DIMS,1)];
    
    V_tilde_tt = E_V{trial_idx}(v_t_t_idx,:);
    P_tilde_tt = E_P{trial_idx}(p_t_t_idx,:);
    P_tilde_ttp1 = P_tilde_tt + V_tilde_tt*DT;
    
    valid_t = ~(any(isnan(P_t),1) | any(isnan(P_tp1),1) | ...
        any(isnan(P_tilde_tt),1) | any(isnan(P_tilde_ttp1),1));
    valid_t(1:T_START-1) = false;
    
    errors.model(trial_idx) = compute_average_absolute_angle(P_tilde_tt(:,valid_t), P_tilde_ttp1(:,valid_t), G, acceptance_zone_radius);
    errors.cursor(trial_idx) = compute_average_absolute_angle(P_t(:,valid_t), P_tp1(:,valid_t), G, acceptance_zone_radius);
end
end


function LL = compute_LL(data, estParams, T_START, DATA_ARE_TRAINING_DATA, MAX_ITERS, VERBOSE)

if VERBOSE
    fprintf('\tEvaluating log-likelihood.\n');
end

if DATA_ARE_TRAINING_DATA
    LL = velime_LL(data, estParams, 'T_START', T_START);
else
    % Evaluating the log-likelihood of the test data requires finding
    % test-data-specific alpha(t). To acquire these, run EM with fixed
    % values of the internal model parameters (A,B,b0,W_v}.
    
    % First remove the training-data-specific alpha (and corresponding
    % noise variance, R). Otherwise, EM will try to apply the training
    % alpha to the test data, which will throw if the training and test
    % datasets are different sizes.
    init_testParams = rmfield(estParams,{'alpha','R'});
    
    % Run EM, only learning the alpha(t).
    TAU = estParams.TAU;
    testParams = velime_fit(data, TAU, ...
        'INIT_METHOD','init_params', ...
        'INIT_PARAMS',init_testParams, ...
        'DO_LEARN_M_PARAMS',false, ...
        'MAX_ITERS', MAX_ITERS, ...
        'VERBOSE', VERBOSE);
    
    % Compute the log-likelihood of the test data using the
    % {A,B,b0,W_v} that were fit to the training data and the {alpha,R}
    % that were fit to the test data. Note: this will differ from
    % test_LL(end) if T_START > (TAU+2), e.g., if sweeping TAU and
    % appropriately evaluating each sweep on exactly the same timesteps
    % in the test data.
    LL = velime_LL(data, testParams, 'T_START', T_START);
end
end

function error_angle = compute_average_absolute_angle(P_tt, P_ttp1, G, radius)

error_angles_from_perimeter = angular_error_from_perimeter(P_tt, P_ttp1, G, radius);
error_angle = mean(abs(error_angles_from_perimeter));

end