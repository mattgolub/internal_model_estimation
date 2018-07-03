% This script demonstrates techniques for model selection in the internal
% model estimation (IME) framework, as described in "Internal models for
% interpreting neural population activity during sensorimotor control," by
% Matthew D. Golub, Byron M. Yu, and Steven M. Chase, eLife 2015. The
% article can be found at https://elifesciences.org/articles/10015/
%
% Please consult the README before running this script.
%
% Here, model selection is performed to probe various settings of TAU, the 
% IME hyperparameter that sets the visual feedback delay in an IME model.
%
% Note that these model selection techniques were not employed in 
% Golub et al., eLife 2015. Rather, in that work we simply fixed the visual
% feedback delay at TAU = 3 timesteps (100 ms) for subject A, and TAU = 4
% (133 ms) for subject C. These settings were chosen based on measurements
% of the subjects' visuomotor latencies (see Figure 2A) that did not rely
% on the IME framework. These measurements were consistent with results
% from model selection in a preliminary verion of the IME framework (see
% Golub et al., ICML, 2013).
%
% A known issue with model selection in the version of IME included in this
% codepack (and used throughout Golub et al., eLife 2015) is that TAU
% affects both the feedback delay and the velocity smoothing in an IME
% model. Concretely, in IME the subject's internal estimates of the current  
% cursor state (position and velocity) are formed based on cursor positions 
% only up until TAU timesteps ago. Intuitively, the value of TAU that 
% matches the subject's true visual feedback delay should yield the best
% model fit, since the subject actually used cursor feedback information up 
% through TAU timesteps ago. However, IME also uses TAU timesteps worth of
% velocity smoothing: current state estimates incorporate TAU timesteps of
% spike counts. Spike counts are inherently noisy and one effect of
% increasing TAU is additional denoising via temporal smoothing. Hence, 
% traditional model selection techniques (as demonstrated below) will tend
% to overestimate TAU, and the magnitude of this bias increases with
% increased noise in the spike count data.
%
% Running this script will, for a representative dataset, fit a an IME
% model for each setting of TAU in {0,1,...,9}. Each IME model is evaluated
% and curves of several assessment metrics are plotted as a function of 
% TAU. See comments in plot_sweep_TAU for further detail.
%
% See also README.md, velime_sweep_TAU, plot_sweep_TAU, example1, example2.
%
% @ Matt Golub, 2018.

clearvars
close all

% Set up the paths to the files included in the codepack.
addpath(genpath('./'))

%% This loads 'data', a data struct from one of the closed-loop BMI
% experiments analyzed in Golub et al, 2015 (dataset A010509). See example1
% for a description of the data.

load('example_data.mat');

%% Options for fitting procedure
CANDIDATE_TAUs = 0:9; % measured in timesteps
MAX_ITERS = 100;
DO_CROSS_VALIDATE = false;
VERBOSE = true;

%% Fit IME for each candidate setting of TAU
sweep = velime_sweep_TAU(data,CANDIDATE_TAUs,...
    'DO_CROSS_VALIDATE', DO_CROSS_VALIDATE, ...
    'MAX_ITERS', MAX_ITERS, ...
    'VERBOSE', VERBOSE);

%% Plot model assessments as a function of TAU. The assessments shown are:
% 1) IME log-likelihood
% 2) Angular error
% 3) The average diagonal element of A_v in the extracted internal models.
% See comments in plot_sweep_TAU for further detail.
plot_sweep_TAU(sweep);