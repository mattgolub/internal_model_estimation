function [SPIKES, CURSOR, Xtarget] = extractData_model4(data, varargin)
%   INPUTS:
%       data is a struct array.
%
% OUTPUTS:
% For each trial, t=1:T,
% SPIKES{t}    - spike count matrix (uDim x N_t)
% CURSOR{t} - cursor state matrix (xDim x N_t)
% Xtarget{t} - target position for trial t (2x1)
%
% *****************************************
% ************* TIME INDEXING *************
% *****************************************
% CURSOR{t} adheres to LDS indexing, ie with state = [position; velocity]
% position(t+1) = position(t) + velocity(t)*dt
%
% SPIKES and CURSOR_vel are NOT aligned
% let p = CURSOR{i}(1:xPosDim,:)
% let u = SPIKES{i}
% if SPIKES are smoothed, then p(:,t+1) = p(:,t) + M2*u(t+1) + m0

SMOOTH_SPIKES = false; % boxcar smooth spike counts
USE_VELOCITY = false; % set state = [pos_x; pos_y; vel_x; vel_y]
USE_ACCELERATION = false; % needed for Batista lab data
RECENTER = false; % use target-centric reference frame
BOXCAR_LEN = 5; % only used if SMOOTH_SPIKES==true
USE_FAILED_TRIALS = false;
USE_PVA_QUIRK = true; % SEE COMMENT IN smoothing function below
assignopts(who,varargin);

num_trials = numel(data);

% cell2mat needs cell array to be 1 x num_trials
CURSOR = cell(1,num_trials);
SPIKES = cell(1,num_trials);
Xtarget = cell(1,num_trials);

success = arrayfun(@(d)(d.success),data)';

for trial = 1:num_trials
    if USE_ACCELERATION
        CURSOR{1,trial} = [data(trial).cursor_pos';
            data(trial).cursor_vel';
            data(trial).cursor_acc'];
    elseif USE_VELOCITY
        xVelDim = size(data(1).cursor_vel',1);
        CURSOR{1,trial} = [data(trial).cursor_pos';
            data(trial).cursor_vel'];
    else
        CURSOR{1,trial} = data(trial).cursor_pos';
    end
    
    Xtarget{1,trial} = data(trial).target';
    
    if RECENTER
        % Set reference frame such that pos = (0,0) is the goal for each
        % trial
        CURSOR{1,trial} = bsxfun(@minus, CURSOR{1,trial}, Xtarget{1,trial});
        Xtarget{1,trial}(:) = 0;
    end
    
    if SMOOTH_SPIKES
        if USE_PVA_QUIRK
            % SPIKES{1,trial} = boxcar_smooth_pva_quirk(data(trial).spike_counts,BOXCAR_LEN)';
            SPIKES{1,trial} = data(trial).smoothed_spike_counts';
        else
            SPIKES{1,trial} = boxcar_smooth_rescale_start(data(trial).spike_counts,BOXCAR_LEN)';
        end
    else
        SPIKES{1,trial} = data(trial).spike_counts';
    end
end

if ~USE_FAILED_TRIALS
    SPIKES = SPIKES(success);
    CURSOR = CURSOR(success);
    Xtarget = Xtarget(success);
end

end

function smooth_data = boxcar_smooth_pva_quirk(data, filtlen)
% Smooths down columns in the case of matrix data

% CAREFUL TO BE SURE THAT THE SPIKES BEING PASSED IN ARE FROM THE VERY
% BEGINNING OF THE TRIAL, OTHERWISE YOU ARE APPLYING THIS ZERO-PAD QUIRK TO
% THE WRONG SEQUENCE OF SPIKES, AND YOUR DECODE WILL NOT AGREE WITH THE
% ONLINE DECODE

smooth_data = zeros(size(data));

if ~isempty(data)
    rescale = filtlen./(1:filtlen-1)';
    for i = 1:size(data,2)
        smooth_data(:,i) = conv(data(:,i),[zeros(filtlen,1); ones(filtlen,1)/filtlen], 'same');
    end
end
end

function smooth_data = boxcar_smooth_rescale_start(data, filtlen)
% Smooths down columns in the case of matrix data

smooth_data = zeros(size(data));

if ~isempty(data)
    rescale = filtlen./(1:filtlen-1)';
    for i = 1:size(data,2)
        smooth_data(:,i) = conv(data(:,i),[zeros(filtlen,1); ones(filtlen,1)/filtlen], 'same');
        smooth_data(1:filtlen-1,i) = smooth_data(1:filtlen-1,i).*rescale;
    end
end
end