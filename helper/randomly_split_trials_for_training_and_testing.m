function [train, test] = randomly_split_trials_for_training_and_testing(U,P,TARGETS,n_test_trials_per_target)

num_trials = numel(TARGETS);
[unique_targets,~,trial_target_idx] = unique([TARGETS{:}]','rows');
num_unique_targets = size(unique_targets,1);

test_trial_idx = nan(num_unique_targets, n_test_trials_per_target);
for this_target_idx = 1:num_unique_targets
    % Generate list of trials for this target
    target_trial_numbers = find(trial_target_idx==this_target_idx);
    
    % Count and make sure there are enough trials to this target
    num_trials_this_target = length(target_trial_numbers);
    if num_trials_this_target<n_test_trials_per_target
        this_target = unique_targets(this_target_idx,:);
        error('%d exceeds the trial count (%d) for the target at (%f,%f)',...
            n_test_trials_per_target, num_trials_this_target, this_target(1), this_target(2));
    end
    
    % Shuffle the list of trials for this target
    shuffled_target_trial_numbers = target_trial_numbers(randperm(num_trials_this_target));
    
    % Randomly select test trials for this target
    test_trial_idx(this_target_idx,:) = shuffled_target_trial_numbers(1:n_test_trials_per_target);
end
test_trial_idx = test_trial_idx(:);
train_trial_idx = setdiff(1:num_trials,test_trial_idx);

train.P = P(train_trial_idx);
train.U = U(train_trial_idx);
train.TARGETS = TARGETS(train_trial_idx);

test.P = P(test_trial_idx);
test.U = U(test_trial_idx);
test.TARGETS = TARGETS(test_trial_idx);