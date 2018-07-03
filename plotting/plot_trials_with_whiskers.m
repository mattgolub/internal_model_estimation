function plot_trials_with_whiskers(P, TARGETS, E_P, E_V, estParams, CURSOR_RADIUS, TARGET_RADIUS)
% plot_trials_with_whiskers(P, TARGETS, E_P, E_V, estParams, CURSOR_RADIUS, TARGET_RADIUS)
%
% @ Matt Golub, 2018

clf; hold on;

num_test_trials = numel(P);
for n = 1:num_test_trials
    % Plot target with enlarged radius to match the "acceptance zone" radius
    % from this experiment
    fill_circle(TARGETS{n},TARGET_RADIUS,'g');
    plot_circle(TARGETS{n},CURSOR_RADIUS + TARGET_RADIUS,'color','g');
    
    % Plot the actual cursor trajectory in black
    plot(P{n}(1,:),P{n}(2,:),'k-o');
    
    % At each timestep, plot a "whisker" representing the evolution of the
    % subject's internal estimates of cursor position, starting at the most
    % recently available visual feedback and up through the predicted cursor
    % position at the upcoming timestep.
    T = size(P{n},2);
    for t = (TAU+1):T
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