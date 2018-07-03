function plot_angular_aiming_errors(angular_errors)
% plot_angular_aiming_errors(angular_errors)
%
% INPUTS:
%
%   angular_errors is a struct produced by velime_evaluate and contains the
%   following fields:
%       .cursor: [1 x # trials] vector containing within-trial averages of
%                absolute angular errors in the actual cursor trajectory.
%       .model:  [1 x # trials] vector containing within-trial averages of 
%                absolute angular errors according to an extracted internal
%                model.
%
% see also velime_evaluate
%
% @ Matt Golub, 2018.

xx = [1 2];

[cursor_avg_error, cursor_sem_error] = meanAndSEM(angular_errors.cursor);
[model_avg_error, model_sem_error] = meanAndSEM(angular_errors.model);
plot_bar_and_sem(xx,[cursor_avg_error model_avg_error],[cursor_sem_error model_sem_error]);

xlim(xx + 0.5*[-1 1]);
ylabel('absolute angular error (degrees)');
set(gca,'xtick',xx,'xticklabel',{'cursor','internal model'},'tickdir','out');
box off

end

function [m, sem, n] = meanAndSEM(x)

x = x(:);
m = nanmean(x);
n = sum(~isnan(x));
sem = nanstd(x,1)/sqrt(n);

end

function B = plot_bar_and_sem(x,y,sem)

BAR_WIDTH = 0.8;

B = bar(x,y,BAR_WIDTH);
hold on

for i = 1:length(x)
    line(x([i i]),[y(i)-sem(i) y(i)+sem(i)],'color','k');
end

set(B,'facecolor','none');

end