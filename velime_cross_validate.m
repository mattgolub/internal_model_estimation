function [estParams, predictions, evaluations, cv_folds] = velime_cross_validate(data, TAU, varargin)
% [estParams, predictions, evaluations, cv_folds] = velime_cross_validate(data, TAU, ...)
%
% Randomly splits data into K blocks, where each block contains no more 
% than one trial to each target. Then for each block, create a disjoint
% training and testing sets of trials. For a given block, the training
% trials contain all trials except that block, and the testing trials
% contain only that block. Then, for each paired set of training and
% testing trials:
%   --extract an internal model (estParams) fit to the training trials
%   --generates cross-validated predictions of the subject's internal state
%     estimates in the testing trials.
%   --computes the angular errors of those cross-validated predictions
%     along with corresponding angular errors in the actual cursor
%     trajectories.
% 
% See Golub, et al., eLife, 2015: "Materials and methods: Computing cross-
% validated internal model predictions"
%
% INPUTS:
%
% data: 
%   Struct with the following fields:
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
% TAU: 
%   The sensory feedback delay. [scalar non-negative integer]. Measured in
%   timesteps.
%
% OUTPUTS:
%
% estParams:
%   K-element struct array, where each element contains the internal
%   model parameters fit to a different set of training trials. Each struct
%   element is as described in velime_fit.
%
% predictions:
%   Struct containing [1 x # trials] cell arrays E_P and E_V, which contain
%   cross-validated predictions of the subject's internal position (E_P)
%   and velocity (E_V) estimates. E_P{i} and E_V{i} contain predictions for
%   trial i in data (i.e., corresponding to data.spike_counts{i}, 
%   data.cursor_position{i}, and data.target_position{i}). Predictions for
%   trial i are generated based on an internal model that was fit to data 
%   not containing trial i. Each element of E_P and E_V are as described in
%   detail in velime_predict.
%
% evaluations:
%   A struct containing the following fields:
%   
%   angular_error.model and angular_error.cursor:  
%       [1 x # trials] vectors containing cross-validated angular errors 
%       according to the subject's internal model (.model) and from the 
%       actual cursor trajectory (.cursor). Element i of each vector 
%       represents an average error across timesteps within trial i. 
%       Internal model errors are based on the cross-validated predictions
%       described above. Additional detail is provided in velime_evaluate.   
%
%   LL: 
%       Scalar value indicating the cross-validated log-likelihood of the
%       data. Value is only computed if optional input DO_COMPUTE_LL is set
%       to true. See velime_evaluate for additional detail.
%
% cv_folds:
%   K-element cell array indicating the blocks of trials used in the cross-
%   validation procedure. Element k is an array containing indices
%   corresponding to a testing set of trials for which predictions and 
%   evaluations were generated based on the internal model in estParams(k)
%   (which was fit to all trials except those indexed by cv_folds{k}). Each
%   trial in data is indexed exactly once across all arrays in cv_folds.
%   See generate_shuffled_blocks for more additional details.
%
% OPTIONAL INPUTS (and default values):
%
%   VERBOSE (false):
%       Logical indicating whether or not to print status updates to the
%       screen.
%
%   T_START (TAU + 2):
%       Integer indicating the timestep to begin evaluating the held-out
%       test data (see velime_evaluate).
%
%   DO_COMPUTE_LL (false):
%       Logical indicating whether or not to compute the cross-validated
%       log-likelihood of the data (see velime_evaluate). If true,
%       substantial additional computation is introduced to this procedure.
%
%   Additional optional arguments specific to velime_fit that are passed in
%   will be passed along to velime_fit (e.g., MAX_ITERS).
%
% See also velime_fit, velime_predict, velime_evaluate, example2, example3,
% generate_shuffled_blocks.
%
% @ Matt Golub, 2018.

VERBOSE = false;
T_START = TAU + 2;
DO_COMPUTE_LL = false;
velime_fit_args = assignopts(who,varargin);

cv_folds = generate_shuffled_blocks(data.target_position);
num_cv_folds = numel(cv_folds);
num_trials = numel(data.cursor_position);

% To store cross-validated internal state estimates
E_P = cell(1,num_trials); E_V = cell(1,num_trials);

% To store cross-validated angular errors
angular_error.model = nan(1,num_trials);
angular_error.cursor = nan(1,num_trials);

% To store cross-validated log-likelihoods
LL_fold = nan(1,num_cv_folds);

for fold_idx = 1:num_cv_folds
    if VERBOSE
        fprintf('Beginning cross-validation fold %d of %d.\n',fold_idx,num_cv_folds);
    end
    
    test_idx = cv_folds{fold_idx};
    test_data = subsample_trials(data, test_idx);
    
    all_folds_but_one = setdiff(1:num_cv_folds,fold_idx);
    train_idx = cell2mat(cv_folds(all_folds_but_one));
    train_data = subsample_trials(data, train_idx);
    
    % Fit velocity-IME model
    estParams(fold_idx) = velime_fit(train_data, TAU, ...
        'VERBOSE', VERBOSE, ...
        velime_fit_args{:});
    
    if VERBOSE
        fprintf('\tExtracting cross-validated predictions.\n');
    end
    
    % Extract prior latent variable distributions ("whiskers")
    [E_P(test_idx), E_V(test_idx)] = velime_predict(test_data, estParams(fold_idx));
    
    if VERBOSE
        fprintf('\tEvaluating cross-validated predictions.\n');
    end
    
    [fold_angular_errors, LL_fold(fold_idx)] = velime_evaluate(test_data, ...
        E_P(test_idx), E_V(test_idx), estParams(fold_idx), ...
        'T_START', T_START, ...
        'DO_COMPUTE_LL', DO_COMPUTE_LL, ...
        'VERBOSE', VERBOSE);
    
    angular_error.model(test_idx) = fold_angular_errors.model;
    angular_error.cursor(test_idx) = fold_angular_errors.cursor;
        
    if VERBOSE
        fprintf('Done.\n\n');
    end
end

predictions.E_P = E_P;
predictions.E_V = E_V;
evaluations.angular_error = angular_error;
evaluations.LL = sum(LL_fold);