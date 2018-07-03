function [M1, M2, m0, V, mdpt] = current_regression(data)
% [M1, M2, m0, V, mdpt] = current_regression(data)
%
% For position only state: initialize M2, m0 by assuming M1 = I and assume
% the subject aims from the current cursor positino straight to the target
% with a constant speed corresponding to the mean distance per timestep.
%
% @ Matt Golub, 2014.

U = data.spike_counts;
X = data.cursor_position;
Xtarget = data.target_position;

POS_IDX = 1:2;

num_trials = numel(U);

mdpt = mean_distance_per_timestep(X);

straight2target = cell(1,num_trials);
straight2target_direction = cell(1,num_trials);
intended_velocity = cell(1,num_trials);
for trialNo = 1:num_trials
    % Unit vectors pointing from Xt to target
    straight2target{trialNo} = bsxfun(@minus,Xtarget{trialNo},X{trialNo}(POS_IDX,:));
    straight2target_norm = sqrt(sum(straight2target{trialNo}.^2,1));
    straight2target_direction{trialNo} = bsxfun(@times,straight2target{trialNo},1./straight2target_norm);
    
    % Scale unit vectors so lengths are the mean
    % distance traveled per timestep
    intended_velocity{trialNo} = mdpt * straight2target_direction{trialNo};
end

[M2, m0, V] = regress_with_delay(intended_velocity,U,-1); % x_t not yet displayed when u_t recorded
M1 = eye(length(m0));
V = diag(V);