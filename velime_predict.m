function [E_P, E_V] = velime_predict(data, estParams)
% [E_P, E_V] = velime_predict(data, estParams)
%
% Extract the subject's internal state predictions given available visual 
% feedback of cursor position and previously issued neural activity. 
% Importantly, these predictions do not take into account target positions.
% Mathematically, these predictions are described by Equation 15 in the
% paper: the prior expectation of the subject's internal state predictions.
% (posterior distributions, as used during model fitting but not here, 
% take into account target positions).
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
% estParams: 
%   struct containing IME parameters as identified by velime_fit. See 
%   velime_fit for details.
%
% OUTPUTS:
%
% E_P and E_V are each [1 x # trials] cell arrays. Each element is
% [2*(TAU+1) x # timesteps]. The details that follow use the notation:
% P = data.cursor_position; U = data.spike_counts.
%
% E_P{n}(:,t) is a column composed of all 2D internal position esimates
% made by the subject at timestep t in trial n. The first of these,
% E_P{n}(1:2,t), is the actual cursor position P{n}(:,t-TAU). Together,
% these positions define a single "whisker". Mathematically, E_P{n}(:,t)
% corresponds to:
%
%    [ \tilde{p}_{t-tau}^t ]                   [ u_{t-tau+1} ]
% E( |        ...          | given p_{t-tau} , |    ...      | )
%    [   \tilde{p}_t^t     ]                   [    u_t      ]
%
% where \tilde{p}_k^t is the subject's internal estimate of the cursor
% position at timestep k given the neural activity up to and including
% timestep t (see Equations 9-11 and Figure 3-figure supplement 1).
%
% Similarly, E_V{n}(:,t) is a column composed of all 2D internal velocity
% estimates made by the subject at timestep t in trial n. Mathematically,
% E_V{n}(:,t) corresponds to:
%
%     [ \tilde{v}_{t-tau}^t ]                     [ u_{t-tau+1} ]
% E ( |         ...         | given v_{t-tau-1} , |    ...      | )
%     [   \tilde{v}_{t}^t   ]                     [   u_{t+1}   ]
%
% where v_{t-tau-1} = (p_{t-tau} - p_{t-tau-1}) / \Delta (Equation 5). The
% E_P{n}(3:end,t) result from integrating the E_V{n}(1:end-2,t) from the
% initial condition E_P{n}(1:2,t) = P{n}(:,t-TAU).
%
% Additional technical notes that may be useful when applying IME to your
% different datasets:
%
% 1) The final velocity estimate, E_V{n}(end-1:end,t), represents the
%    subject's intended velocity command at timestep t, which the subject
%    internally predicts would drive the cursor to position
%    \tilde{p}_{t+1}^t (which does not appear in E_P).
%
% 2) Note that the actual cursor velocity v_{t-tau} (proportional to
%    p_{t-tau} - p_{t-tau-1}) is not included in E_V (whereas the actual 
%    cursor position p_{t-tau} is included in E_P).
%
% 3) Some identities for intuition about the relationships between various
%    data (U, P), extracted model parameters (estParams), and internal
%    state predictions (E_P and E_V) (these hold true for all timesteps
%    where relevant data are defined).
%
%    Using A, B, b0, and dt from estParams, trial number n, and timestep t:
%
%    Define actual cursor velocity: V{n} = diff(P{n},[],2)/dt. Rearranging
%    gives P(:,t) = P(:,t-1) + V(:,t-1) * dt, which matches the convention
%    defined in Equation 5.
%    
%    Internal estimates of cursor velocity:
%       E_V{n}(1:2,t) = A*V{n}(:,t-TAU-1) + B*U{n}(:,t-TAU+1) + b0
%       E_V{n}(3:4,t) = A*E_V{n}(1:2,t) + B*U{n}(:,t-TAU+2) + b0
%                   ...
%       E_V{n}(end-1:end,t) = A*E_V{n}(end-3:end-2,t) + B*U{n}(:,t+1) + b0
%
%    Internal estimates of cursor position:
%       E_P{n}(1:2,t) = P{n}(:,t-TAU)
%       E_P{n}(3:4,t) = E_P{n}(1:2,t) + E_V{n}(1:2,t) * dt
%                   ...
%       E_P{n}(end-1:end,t) = E_P{n}(end-3:end-2,t) + E_V{n}(end-3:end-2,t) * dt
%
% 4) The identites above show that the temporal alignment conventions of the
%    input data are consistently maintained in the model predictions.
%    Concretely, in the data, p(t) = P{n}(:,t) is the first actual cursor
%    position to reflect the spike counts u(t) = U{n}(:,t). Similarly, at
%    timestep t, \tilde{p}_t^t is the first internal prediction to take
%    into account u(t), and similar for \tilde{p}_{t-k}^t and corresponding
%    u(t-k). As IME "unrolls" a new chain of internal state estimates at
%    each timestep, u(t) becomes the most recent spike count to be
%    incorporated into \tilde{p}_t^t (as described above),
%    \tilde{p}_{t+1}^{t-1}
%
% Technical note: this function merely translates predictions from the 
% notational conventions of the fitting code (velime_fit) into the more
% readable notational conventions of the paper. This function is not called
% during fitting (within velime_fit).
%
% See also velime_fit, velime_evaluate, plot_trials_with_whiskers.
%
% @ Matt Golub, 2018

TAU = estParams.TAU;
dt = estParams.dt;

% Get predictions from the model (using the notational convention of the
% fitting implementation)
[C,G,const] = velime_assemble_data(data,TAU,dt);
E_X = velime_prior_expectation(C, estParams);

% Tranlate predictions into the notational conventions of the paper
N_trials = numel(const.trial_map);
pIdx = const.x_idx(1:const.pDim,:);
vIdx = const.x_idx(const.pDim+1:const.xDim,:);
for trialNo = 1:N_trials
    trial_indices = const.trial_map{trialNo};
    T_data = size(data.cursor_position{trialNo},2);
    
    % Pad with TAU + 1 columns of nan's at the start and one column at the
    % end. This ensures that these quantities align in time with the
    % conventions in the input P.
    if T_data >= TAU + 2
        E_P{trialNo} = [nan(numel(pIdx),TAU+1) E_X(pIdx,trial_indices) nan(numel(pIdx),1)];
        E_V{trialNo} = [nan(numel(vIdx),TAU+1) E_X(vIdx,trial_indices) nan(numel(vIdx),1)];
    else
        % Make sure nan padding does not change the length of the trial.
        E_P{trialNo} = nan(numel(pIdx),T_data);
        E_V{trialNo} = nan(numel(vIdx),T_data);
    end
end