function plot_sweep_TAU(sweep)
% plot_sweep_TAU(sweep)
%
% Plot various assessments for IME models fit across a range of settings of
% TAU. Curves are plotted for the following model assessments:
%
% 1) IME log-likelihood (of the target positions given available cursor
%    position feedback and causal neural activity). Larger values indicate
%    better goodness-of-fit.
% 2) Angular errors according the extracted internal models. Smaller values
%    indicate better goodness-of-fit.
% 3) The average diagonal element of A_v in the extracted internal models.
%    A_v controls how the subject's internal state estimates are influenced
%    by velocity feedback (from TAU timesteps ago) and by velocity 
%    smoothing across timesteps. "Larger" values of A_v (here assessed via 
%    the average of the diagonal elements) correspond to IME models that
%    rely more heavily on velocity feedback, which might be an indication 
%    that TAU matches the subject's true visuomotor latency. Note, this
%    assessment should be interpreted with care because A_v also controls
%    the amount of velocity smoothing leveraged by an IME model. More
%    smoothing can be necessary when spike count data are noisy (e.g., due
%    to limited numbers of recorded neurons).
%
% INPUTS:
% 
% sweep: 
%    A struct array produced by velime_sweep_TAU. Each element of sweep
%    corresponds to a different setting of TAU. 
%
% See also velime_sweep_TAU, example3.
%
% @ Matt Golub, 2018

TAUs = [sweep.TAU];

% Find optimal TAU according to...

% log-likelihood
[max_LL, idx_max_LL] = max([sweep.LL]);
TAU_star.LL = TAUs(idx_max_LL);

% angular error
angular_error = arrayfun(@(s)(nanmean(s.angular_error.model)),sweep);
[min_angular_error, idx_min_angular_error] = min(angular_error);
TAU_star.angular_error = TAUs(idx_min_angular_error);

% the autoregressive component of velocity predictions
if numel(sweep(1).estParams) == 1
    % Handles format of estParams that were fit WITHOUT cross validation
    ar_coeffs = arrayfun(@(s)(mean(diag(s.estParams.A))),sweep);
else
    % Handles format of estParams that were fit WITH cross validation
    ar_coeffs = nan(1,numel(sweep));
    for sweep_idx = 1:numel(sweep)
        estParams_array = sweep(sweep_idx).estParams;
        ar_coeffs_per_fold = arrayfun(@(estParams)(mean(diag(estParams.A))),estParams_array);
        ar_coeffs(sweep_idx) = mean(ar_coeffs_per_fold);
    end
end
[max_ar_coeff, idx_max_ar_coeff] = max(ar_coeffs);
TAU_star.ar_coeff = TAUs(idx_max_ar_coeff);

subplot(3,1,1); hold on;
plot(TAUs,[sweep.LL])
plot(TAU_star.LL,max_LL,'ro');
ylabel('Log-likelihood');
box off
set(gca,'tickdir','out')

subplot(3,1,2); hold on;
plot(TAUs,angular_error)
plot(TAU_star.angular_error,min_angular_error,'ro');
ylabel({'Absolute angular','error (degrees)'});
xlabel('\tau');
box off
set(gca,'tickdir','out')

subplot(3,1,3); hold on;
plot(TAUs,ar_coeffs)
plot(TAU_star.ar_coeff,max_ar_coeff,'ro');
ylabel({'Approx autoregressive','coefficient [mean(diag((A_v))]'});
xlabel('\tau');
box off
set(gca,'tickdir','out')

end