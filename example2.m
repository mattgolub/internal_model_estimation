% This script demonstrates a proper cross-validated application of 
% internal model estimation (IME) as described in "Internal models for 
% interpreting neural population activity during sensorimotor control," by 
% Matthew D. Golub, Byron M. Yu, and Steven M. Chase, eLife 2015. The 
% article can be found at https://elifesciences.org/articles/10015/
% 
% Please consult the README before running this script.
% 
% Running this script will, for a representative dataset, fit a set of
% cross-validated IME models, extract the subject's (cross-validated)
% internal estimates of cursor position and velocity, and evaluate errors 
% in the actual cursor trajectory and in the subject's internal estimates. 
% 
% This script generates two figures: one quantifying the cross-validated 
% angular errors (in the style of Figure 3C); and a second showing an 
% example cursor trajectory overlaid with the subject's (cross-validated) 
% internal state estimates (in the style of Figure 4A,B).
%
% See also README.md, velime_cross_validate, example1, example3.
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

%% Options for fitting procedure ------------------------------------------

% Set the random seed so cross validation results are reproducible
% (the only stochastic operation is the partitioning of data into training
% and testing folds).
rng(0) 

VERBOSE = true; % Whether or not to print status updates throughout fitting

MAX_ITERS = 100; % Recommendation: set to 5000 (iterations) unless you know
% you can get away with fewer. Here we use fewer iterations so the example 
% runs faster.
%
% Below are approximate runtimes (3.1Ghz Intel core i7) for this script and
% resulting cross-validated angular errors for the example data and the
% fixed seeding of the random number generator:
%                                            CROSS-VALIDATED
% MAX_ITERS     APPROX RUNTIME (s)       ANGULAR ERROR (degrees)
%    10               2.2                        4.2634
%    25               4.5                        2.9826
%    50               8.3                        2.6319
% * 100 * * * * * *  15.9 * * * * * * * * * * *  2.5827 * * * * * * * * * *                       
%   500              76.1                        2.6016                        
%  1000             163.2                        2.6069                        
%  5000             751.8                        2.6384               

TAU = 3; % Sensory feedback delay, measured in timesteps [scalar non-
% negative integer]. Here, as in the paper, we fix TAU = 3 to be consistent
% with model-free estimates of this animal's sensory feedback delay (100ms,
% see Figure 2A, left). One can fit this as a hyperparameter if desired 
% (see example3.m). 

%% Extracting an internal model -------------------------------------------

% Run fitting procedure and generate cross-validated predictions.
% This uses leave-one-block-out cross validation, where a block contains
% one trial to each of 16 targets. Trials are assigned to blocks randomly
% without replacement. The example data contain 11 blocks of trials, so
% this is K-fold cross validation, where K = 11. See the paper section: 
% "Materials and Methods, "Computing cross-validated internal model 
% predictions".

[estParams, predictions, evaluations, cv_folds] = velime_cross_validate(data, TAU,...
    'MAX_ITERS',MAX_ITERS,'VERBOSE',VERBOSE);

%% Figures ----------------------------------------------------------------

% Generate angular error plot in the style of Figure 3C. 
figure(1); clf;
plot_angular_aiming_errors(evaluations.angular_error);

% Generate example "whisker" plot in the style of Figure 4A,B. Each whisker
% represents the subject's internally predicted evolution of the cursor
% trajectory given the most recently available visual feedback of cursor 
% position and the subsequently issued neural commands.
figure(2); clf;
EXAMPLE_IDX = 98;
ex_data = subsample_trials(data, EXAMPLE_IDX);

% Find the test fold containing the example trial. Then use the parameters
% that were fit to the data in all folds except that held-out test fold.
ex_cv_fold_idx = find(cellfun(@(cv)(any(cv==EXAMPLE_IDX)),cv_folds)); 
plot_trials_with_whiskers(ex_data, ...
    predictions.E_P(EXAMPLE_IDX), predictions.E_V(EXAMPLE_IDX), ...
    estParams(ex_cv_fold_idx));