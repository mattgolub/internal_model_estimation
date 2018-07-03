function [M2, m0, V] = target_regression(data)
% [M2, m0, V] = target_regression(data)
%
% Initialize M2, m0 by fitting
% u_t = M2*d_t + m0 + v_t, v_t ~ N(0,V),
% where d_t is a unit vector pointing from the center of the workspace to 
% the target, scaled by the mean distance traveled per timestep.
%
% Assume desired change in position (i.e. velocity) is given by target
% direction.  Assume desired change in higher order kinematics is zero
% (i.e. acceleration, jerk).
%
% @ Matt Golub, 2014.

U = data.spike_counts;
X = data.cursor_position;
Xtarget = data.target_position;

xDim = size(X{1},1);
num_trials = numel(U);
T = cellfun(@(u)(size(u,2)),U);

distances = cellfun(@(x)(sqrt(sum(diff(x,1,2).^2,1))),X,'uniformoutput',false);
mean_distance_per_timestep = mean(cell2mat(distances));

target_direction = cell(1,num_trials);
intended_velocity = cell(1,num_trials);
intended_kinematics = cell(1,num_trials);
for trialNo = 1:num_trials
    % Unit vectors pointing from center to target
    target_direction{trialNo} = repmat(Xtarget{trialNo},1,T(trialNo));
    target_norm = sqrt(sum(target_direction{trialNo}.^2,1));
    target_direction{trialNo} = bsxfun(@times,target_direction{trialNo},1./target_norm);
    
    % Scale unit vectors so lengths are the mean
    % distance traveled per timestep
    intended_velocity{trialNo} = mean_distance_per_timestep * target_direction{trialNo};
    
    intended_kinematics{trialNo} = [intended_velocity{trialNo}; zeros(xDim-2,size(intended_velocity{trialNo},2))];
end

[M2, m0, V] = regress_with_delay(intended_kinematics,U,0);
V = diag(V);