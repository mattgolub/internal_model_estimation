function H = fill_circle(center, radius, color, varargin)
% H = fill_circle(center, radius, color, varargin)
%
% center is 1x2, radius is 1x1
% H is a handle to the plotted circle
%
% @ Matt Golub, 2014.

center = reshape(center,1,2);

N = 3000;
r = linspace(0,2*pi,N);
unit_circle = [cos(r)' sin(r)'];
target_pts = bsxfun(@plus, radius*unit_circle, center);

% Line target
if ~ishold(gca)
    cla
end
H = fill(target_pts(:,1), target_pts(:,2),...
    color,'linestyle','none',varargin{:});