% This script runs a minimal working example of internal model estimation
% (IME) as described in "Internal models for interpreting neural population
% activity during sensorimotor control," by Matthew D. Golub, Byron M. Yu, 
% and Steven M. Chase, eLife 2015. The article can be found at
% https://elifesciences.org/articles/10015/
% 
% Please consult the README before running this script.
% 
% Running this script will, for a representative dataset, fit a single IME 
% model, extract the subject's internal estimates of cursor position and 
% velocity, and evaluate errors in the actual cursor trajectory and in the 
% subject's internal estimates. This script will also generate two figures:
% one quantifying those errors; and a second showing a cursor trajectory 
% and the subject's internal state estimates from an example trial.
%
% IMPORTANT! WHAT FOLLOWS IS NOT AN APPROVED USAGE OF IME. This example 
% generates predictions using the same data used to extract the internal 
% model. We include this minimal working example to illustrate the key 
% components of this codepack in a script that runs very quickly. 
% As with any machine learning procudure, meaningful interpretation of
% model predictions requires cross validation (i.e., generating predictions
% using "test" or "validation" data that was held out during model fitting.
% We demonstrate an example of such an approved usage of IME in example2.m.
%
% See also README.md, velime_fit, velime_predict, velime_evaluate, 
% example2, example3.
% 
% @ Matt Golub, 2018.

clearvars
close all

% Set up the paths to the files included in the codepack.
addpath(genpath('./'))

%% This loads 'data', a data struct from one of the closed-loop BMI 
% experiments analyzed in Golub et al, 2015 (dataset A010509).

load('example_data.mat');

% The struct contains the following fields:
%
% spike_counts: Binned spike counts. [1 x # trials] cell array. Each 
% element is [# neurons x # timesteps]. In these data, each timestep is
% 33ms, and bins were non-overlapping.
%
% cursor_position: 2D cursor positions. [1 x # trials] cell array. Each 
% element is [2 x # timesteps].
%
% target_position: 2D target positions. [1 x # trials] cell array. Each 
% element is [2x1]. The fitting code will also accept time-varying target
% position, in which case this can be [2 x # timesteps].
%
% cursor_radius and target_radius: [non-negative scalars]. These are the
% radii of the circular cursor and targets, respectively. Trial success
% required a visual overlap of the cursor and the target. These values are
% used for computing angular errors (defined relative to the target
% perimeter) and for plotting example trials.
%
% Values in cursor_position, target_position, cursor_radius and 
% target_radius must all be given in the same units. In the example data, 
% those units are mm. The 
%
% The columns of each element of spike_counts and cursor_pos should be 
% aligned such that:
%
%   spike_counts{n} = [u(1) u(2) ... u(T)]
%   cursor_pos{n}   = [p(1) p(2) ... p(T)]
%
% where n indexes trials, p(t) is the first cursor position to take into
% account the spike count vector u(t). For example, this convention is
% satisfied by the following linear decoder update equations:
%
% Example 1: A population vector algorithm or ordinary linear estimator.
%       p(t) = p(t-1) + v(t)*dt
%       v(t) = B*u(t) + b0
%
% Example 2: A velocity Kalman filter:
%       p(t) = p(t-1) + v(t)*dt
%       v(t) = M1*v(t-1) + M2*u(t) + m0

%% Options for fitting procedure ------------------------------------------

TAU = 3; % Sensory feedback delay, measured in timesteps [scalar non-
% negative integer]. One can fit this as a hyperparameter if desired, but 
% here we simply fix it.

MAX_ITERS = 100; % Recommendation: set to 5000 (iterations) unless you know
% you can get away with fewer. Here we use fewer iterations so the example 
% runs faster.
%
% Below are approximate runtimes (3.1Ghz Intel core i7) for this script and
% resulting angular errors for the example data.
%                                            
% MAX_ITERS     APPROX RUNTIME (s)       ANGULAR ERROR (degrees)
%    10               0.3                        3.8766
%    25               0.5                        2.5144
%    50               0.8                        2.1383
% * 100 *           * 1.6 *                   *  2.0677 *
%   500               7.2                        2.0774                        
%  1000              14.4                        2.0837
%  5000              58.0                        2.0991
%
% Note that angular error begins to increase ever so slightly after roughly
% 100 iterations. This is because IME does not directly optimize angular 
% errors. Rather, it optimizes a likelihood whose improvements tend to
% correlate highly with improvements in angular error.
%
% Although this situation is somewhat reminiscent of overfitting, we can
% clearly rule that out in this case, because in this example we are 
% evaluating angular errors on the training data. To assess overfitting, 
% one must evaluate on held out data not used for training; see cross 
% validation in example2.m). Rather, because 

VERBOSE = true; % Whether or not to print status updates throughout fitting

% -------------------------------------------------------------------------

%% Run fitting procedure to extract the subject's internal model of the BMI
estParams = velime_fit(data, TAU, 'VERBOSE', VERBOSE, 'MAX_ITERS', MAX_ITERS);

%% Extract the subject's internal estimates of cursor position and velocity
[E_P, E_V] = velime_predict(data, estParams);

%% Evaluate angular errors in actual cursor trajectories and according to
% the subject's internal model
angular_error = velime_evaluate(data, E_P, E_V, estParams);

%% Figures ----------------------------------------------------------------

% Generate angular error plot in the style of Figure 3C. 
figure(1); clf;
plot_angular_aiming_errors(angular_error);

% Generate example "whisker" plot in the style of Figure 4A,B. Each whisker
% represents the subject's internally predicted evolution of the cursor
% trajectory given the most recently available visual feedback of cursor 
% position and the subsequently issued neural commands.
figure(2); clf;
EXAMPLE_IDX = 98;
ex_data = subsample_trials(data, EXAMPLE_IDX);
plot_trials_with_whiskers(ex_data, E_P(EXAMPLE_IDX), E_V(EXAMPLE_IDX), estParams);