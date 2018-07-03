function H = plot_circle(center, radius, varargin)
% H = plot_circle(center, radius, varargin)
% center is 1x2, radius is 1x1
% H is a handle to the plotted circle.
%
% @ Matt Golub, 2014

center = reshape(center,1,2);

unit_circle = [cos(0:.01:2*pi)' sin(0:.01:2*pi)'];
target_pts = bsxfun(@plus, radius*unit_circle, center);

% Line target
if ~ishold(gca)
    cla
end
H = line(target_pts(:,1), target_pts(:,2),...
    varargin{:});