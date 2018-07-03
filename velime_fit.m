function [estParams, LL] = velime_fit(data, TAU, varargin)
% [estParams, LL] = velime_fit(data, TAU,...)
%
% Fits the internal model estimation (IME) framework via expectation
% maximization (EM). The IME framework is described in Golub, Yu & Chase,
% eLife, 2015 (https://elifesciences.org/articles/10015). In the
% descriptions below, we reference specific equations from this paper.
%
% INPUTS:
%
% data.spike_counts: 
%   Binned spike counts. [1 x # trials] cell array. Each element is
%   [# neurons x # timesteps].
%
% data.cursor_position: 
%   2D cursor positions. [1 x # trials] cell array. Each element is
%   [2 x # timesteps].
%
% data.target_position: 
%   2D target positions. [1 x # trials] cell array. Each element is [2x1]
%   or [2 x # timesteps].
%
% TAU: 
%   The sensory feedback delay. [scalar non-negative integer]. Measured in
%   timesteps.
%
% OUTPUTS: 
%
% estParams: a struct containing the following parameters:
%
%   The subject's internal model of the BMI mapping:
%         A: [2 x 2]. This corresponds to \tilde{A}_v in Equation 10.
%         B: [2 x # neurons]. This corresponds to \tilde{B}_v in Equation 10.
%        b0: [2 x 1]. This corresponds to \tilde{b}_v in Equation 10.
%
%   Isotropic noise variances to account for internal state predictions not
%   captured by the internal model:
%       W_v: [2 x 1]. This corresponds to the variance (w) of the noise
%            variables w^t_k in Equation 10.
%       W_p: [2 x 1]. These will be zeros (because position updates are
%            defined deterministically given the subject's internal
%            estimate of velocity).
%
%   Parameters of the observation model (note that these relate only to the
%   training data and are not used when generating predictions on held-out
%   test data).
%     alpha: [1 x # of "whiskers"]. These the distance scale parameters of
%     Equation 14.
%         R: [2x1 double]. Isotropic noise variances to account for
%         internal velocity predictions that do not point straight to the
%         target from the subject's up-to-date position prediction.
%
% The following are not fit to data and are simply carried around to
% facilitate generating predictions after fitting.
%
%       TAU: [1 x 1]. This is the TAU specified above.
%        dt: Refers to \Delta in Equations 5 and 9. Without loss of
%            generality, this codepack fixes dt = 1. There is no need for
%            dt to match the actual timestep in the data, as dt is
%            only used to compute empirical velocities from the position
%            data, and to integrate velocities to generate position
%            estimates. These reciprocal computations cancel all effects
%            of dt.
%
% LL: a vector containing the trajectory of training log-likelihoods. EM
%     guarantees that these values increase with each iteration.
%
% OPTIONAL INPUTS (and default values):
%
%   TOL (1e-8): 
%       Convergence criteria for across-iteration change in training log-
%       likelihood.
%   ZERO_TOL(1e-10): 
%       Criteria for detecting violations of EM-guaranteed increases in
%       training log-likelihood. Useful for debugging when modifying the
%       algorithm implementation.
%   MAX_ITERS (5e3):
%       The maximum number of EM iterations to perform. The algorithm will
%       terminate upon reaching this iteration number unless convergence
%       was detected before performing this many iterations.
%   MIN_ITERS (10):
%       The minimum number of EM iterations to perform. The algorithm will
%       not terminate in fewer than this many iterations, even if
%       convergence was detected before performing this many iterations.
%   VERBOSE (false):
%       Logical indicating whether or not to print status updates to the
%       screen.
%   INIT_METHOD ('current_regression'):
%       String indicating the initialization method. Valid choices are
%       'current_regression' (detailed in the paper), 'target_regression'
%       (see initialize_velime.m for details), or 'init_params' (see
%       below).
%   INIT_PARAMS (NaN):
%       An initial set of parameters for EM (matching the format of
%       estParams, as described under OUTPUTS above). To use these initial
%       parameters, INIT_METHOD must be set to 'init_params'.
%   DO_LEARN_M_PARAMS (true):
%       Logical indicating whether or not to learn the dynamics parameters
%       {A,B,b0,W_v}. If false, values will be fixed as initialized.
%   DO_LEARN_ALPHA_PARAMS (true):
%       Logical indicating whether or not to learn the scale parameters
%       {alpha, R}. If false, values will be fixed as initialized.

startTime = tic;

initialize_velime

if VERBOSE
    fprintf('\tFitting velocity-IME via EM.\n')
end

LLi   = nan;
LL    = nan(1,MAX_ITERS);
iter_times = nan(1,MAX_ITERS);
iters_completed = 0;

while (MAX_ITERS>0)
    iters_completed_start_time = tic;
    
    LLold = LLi;
    
    % =======
    % E-step
    % =======
    [LLi, E_X_posterior, COV_X_posterior]  = velime_estep(C, G, estParams, const);
    
    LL(iters_completed+1) = LLi;
        
    if mod(iters_completed,50)==0 && iters_completed>0
        if VERBOSE
            fprintf('\t\tvelime(TAU=%d) iters: %d, Mean iter time: %1.1es, LL improvement: %.3e\n', TAU, iters_completed,mean(iter_times(1:iters_completed)), LLi-LLold);
        end
    end
    
    % =======
    % Check for convergence
    % =======
    if iters_completed>1 && LLi < LLold && abs(LLold-LLi)>ZERO_TOL
        fprintf('\t\tVIOLATION velime(TAU=%d) iters: %d: %g\n', TAU, iters_completed, LLi-LLold);
        try
            plot(LL); drawnow;
        catch
            % do nothing
        end
    end
    if (iters_completed>MIN_ITERS) && (LLi>LLold) && (LLi - LLold < TOL)
        if VERBOSE
        fprintf('\t\tvelime(TAU=%d) converged after %.2f seconds: %d iters\n',TAU, toc(startTime),iters_completed);
        end
        break; % convergence
    elseif iters_completed>=MAX_ITERS
        if VERBOSE
            fprintf('\t\tvelime(TAU=%d) iter limit reached, aborting after %.2f seconds, LL improvement: %.3e\n', TAU, toc(startTime), LLi-LLold);
        end
        break
    end
    
    % =======
    % M-step
    % =======
    if DO_LEARN_M_PARAMS
        M_Params_fastfmc = velime_mstep_MParams(E_X_posterior, COV_X_posterior, C, const);
        estParams.A = M_Params_fastfmc.M1;
        estParams.B = M_Params_fastfmc.M2;
        estParams.b0 = M_Params_fastfmc.m0;
        estParams.W_v = M_Params_fastfmc.V;
    end
    
    if DO_LEARN_ALPHA_PARAMS
        ALPHA_Params = velime_mstep_alphaParams(G, E_X_posterior, COV_X_posterior, const);
        estParams.alpha = ALPHA_Params.alpha;
        estParams.R = ALPHA_Params.R;
    end
    
    iters_completed = iters_completed + 1;
    iter_times(iters_completed) = toc(iters_completed_start_time);
end

end