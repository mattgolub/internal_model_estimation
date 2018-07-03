function [block_trial_indices, block_targets] = generate_shuffled_blocks(trials_targets)
% [block_trial_indices, block_targets] = generate_shuffled_blocks(trials_targets)
%
% Generates sets of non-overlapping trial indices with matched target
% distributions for use in cross validation. If trials_targets contains the
% same number of trials to each target, each set of trial indices will be
% the same size and each set will index to exactly one trial per target.
%
% INPUTS:
%   trial_targets: [1 x # trials] cell array, where each element is a
%   column vector indicating a target for that trial.
%
% OUTPUTS:
%   block_trial_indices: K-element cell array containing random,
%   non-overlapping sets of trial indices.
%
%   block_targets: K-element cell array containing the target IDs of the
%   trials indexed by block_trial_indices. Target IDs are self consistent
%   in their relationship to the targets in trials_targets, but are 
%   internally defined arbitrarily. These can be used to verify the
%   distribution of targets indexed by each set in block_trial_indices.
%
% @ Matt Golub, 2014-2018.

% Determine unique targets (targets), which target for each trial
% (trials_target_idx), and which trials for each target (trials_idx_per_target)
[unique_targets, ~, trials_target_idx] = unique([trials_targets{:}]','rows');

% Find all trials with each target
num_targets = size(unique_targets,1);
num_target_repetitions = nan(num_targets,1);
for target_idx = 1:num_targets
    trials_idx_per_target{target_idx,1} = find(trials_target_idx==target_idx);
    num_target_repetitions(target_idx) = length(trials_idx_per_target{target_idx});
    
    % Shuffle trial numbers for this target
    shuffle_idx = randperm(num_target_repetitions(target_idx));
    shuffled_trials_idx_per_target{target_idx,1} = trials_idx_per_target{target_idx}(shuffle_idx);
end
max_num_trials_per_target = max(num_target_repetitions);

% Pad shuffled trial indices with nans so all cells contain an array of the 
% same length.
shuffled_trials_idx_per_target_nanpad = cellfun(@(x)([x; nan(max_num_trials_per_target-length(x),1)]),shuffled_trials_idx_per_target,'uniformoutput',false);
trials_idx_mat = cell2mat(shuffled_trials_idx_per_target_nanpad');

% Make cell array where each cell contains a random trial number to each
% target.  If targets are uniformly represented in data, each resulting
% cell element will be the same length (equal to the number of targets).
for block_idx = 1:size(trials_idx_mat,1)
    temp = trials_idx_mat(block_idx,:);
    block_trial_indices{block_idx,1} = temp(~isnan(temp))';
    block_targets{block_idx,1} = trials_target_idx(block_trial_indices{block_idx,1},:);
end

end