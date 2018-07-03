function [angles, theta, phi] = angular_error_from_perimeter(x_t_t, x_tp1_t, target_center, radius)
% [angles, theta, phi] = angular_error_from_perimeter(x_t_t, x_tp1_t, target_center, radius)
%
% INPUTS:
% x_t_t and x_tp1_t are DxN, N = number of time points, D = internal state
% dimensionality.  Only position elements (1:2,:) are looked at.
%
% target_center is 2x1, radius is 1x1
%
% OUTPUTS: 
% angles is Nx1 signed angular errors from the target perimeter, ie the
% minimal angular rotation of the aiming direction such that rotating
% (x_tp1_t - x_t_t) about x_t_t would intersect with a point, p, on the
% target perimeter
%
% OPTIONAL OUTPUTS:
% theta is the angular error between the aiming direction and (target_center - x_t_t)
% phi is the max angle between (target_center-x_t_t) and a perimeter point
%
% @ Matt Golub, 2014

[D,N] = size(x_t_t);
[D2,N2] = size(x_tp1_t);

if N~=N2 || D~=D2 || numel(target_center)~=2 || size(target_center,1)~=2 || numel(radius)~=1
    error('Inputs are not formatted correctly');
end

pos_idx = [1 2];

v1 = bsxfun(@minus, target_center, x_t_t(pos_idx,:));
v2 = x_tp1_t(pos_idx,:) - x_t_t(pos_idx,:);
v1_dot_v2 = dot(v1,v2);
% This is the sign of the zero-padded cross product:
% v1_cross_v2 = cross([v1;zeros(1,N)],[v2;zeros(1,N)]);
error_is_counterclockwise = sign(v1(1,:).*v2(2,:) - v1(2,:).*v2(1,:));
length_v1 = sqrt(sum(v1.^2,1));
length_v2 = sqrt(sum(v2.^2,1));

% Find the angle between v1 and v2
cos_theta = v1_dot_v2./(length_v1.*length_v2);
theta = acosd(cos_theta);
theta(error_is_counterclockwise==-1) = -theta(error_is_counterclockwise==-1);

% Find the largest angle between v1 and v3,
% where v1 is x_t_t to target_center, and
% v3 is x_t_t to a point on the circle perimeter
phi = asind(radius./length_v1) .* sign(theta);
phi(length_v1<radius) = 180;

radius_sq = radius^2;

% theta(n) and phi(n) must be the same sign.
idx_zero = abs(theta)<abs(phi); % aiming direction would bring cursor to target

angles = nan(1,N);
angles(idx_zero) = 0;
angles(~idx_zero) = theta(~idx_zero)-phi(~idx_zero);