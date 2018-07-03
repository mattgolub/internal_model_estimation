function [U, X] = cell2mat_forLaggedRegression(u, x, d)
% Preprocessing step for lagged regression: 
% u(t) = L*x(t-d) + l0 + noise

num_trials = numel(u);

X_cell = cell(1,num_trials);
U_cell = cell(1,num_trials);
if d>=0
    % x lags u 
    for trialno = 1:num_trials
        % Append x and u data such that x lags behind u
        X_cell{trialno} = x{trialno}(:,1:end-d);
        U_cell{trialno} = u{trialno}(:,d+1:end);
    end
else % u lags x
    dpos = abs(d);
    for trialno = 1:num_trials
        % Append x and u data such that x lags behind u
        X_cell{trialno} = x{trialno}(:,dpos+1:end);
        U_cell{trialno} = u{trialno}(:,1:end-dpos);
    end
end
X = cell2mat(X_cell);
U = cell2mat(U_cell);