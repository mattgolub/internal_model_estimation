function sub_data = subsample_trials(data, trial_idx)
% sub_data = subsample_trials(data, trial_idx)
%
% Returns a subsampled dataset containing only the trials specified in
% trial_idx.
%
% @ Matt Golub, 2018.

sub_data.cursor_position = data.cursor_position(trial_idx);
sub_data.spike_counts = data.spike_counts(trial_idx);
sub_data.target_position = data.target_position(trial_idx);

sub_data.cursor_radius = data.cursor_radius;
sub_data.target_radius = data.target_radius;