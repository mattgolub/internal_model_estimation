function plot_trials_with_whiskers(data, E_P, E_V, estParams)
% plot_trials_with_whiskers(data, E_P, E_V, estParams)
%
% Plot cursor trajectories from the trials in data, overlaid by "whiskers"
% representing the subject's internally predicted evolution of the 
% cursor trajectory given the most recently available visual feedback of
% cursor position and the subsequently issued neural commands.
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
%   velime_predict. Each element is [2*(TAU+1) x # timesteps].
%
% See also velime_fit, velime_predict
%
% @ Matt Golub, 2018

P = data.cursor_position;
TARGETS = data.target_position;
target_radius = data.target_radius;
cursor_radius = data.cursor_radius;

clf; hold on;

num_test_trials = numel(P);
for n = 1:num_test_trials
    % Plot target with enlarged radius to match the "acceptance zone" radius
    % from this experiment
    fill_circle(TARGETS{n},target_radius,'g');
    plot_circle(TARGETS{n},cursor_radius + target_radius,'color','g');
    
    % Plot the actual cursor trajectory in black
    plot(P{n}(1,:),P{n}(2,:),'k-o');
    
    % At each timestep, plot a "whisker" representing the evolution of the
    % subject's internal estimates of cursor position, starting at the most
    % recently available visual feedback and up through the predicted cursor
    % position at the upcoming timestep.
    T = size(P{n},2);
    for t = (estParams.TAU+1):T
        % \tilde{p}_k^t for k = t-TAU, ..., t
        plot(E_P{n}(1:2:end,t),E_P{n}(2:2:end,t),'r.-');
        
        % \tilde{p}_t^t, the subject's up-to-date estimate of the current
        % cursor position (at timestep t)
        p_t_t = E_P{n}(end-1:end,t);
        plot(p_t_t(1),p_t_t(2),'ro');
        
        % Compute \tilde{p}_{t+1}^t, the predicted cursor position after the
        % subject issues the intended velocity command, \tilde{v}_t^t
        v_t_t = E_V{n}(end-1:end,t);
        p_tp1_t = p_t_t + v_t_t * estParams.dt;
        plot([p_t_t(1) p_tp1_t(1)],[p_t_t(2) p_tp1_t(2)],'r-');
    end
end

axis image
axis off