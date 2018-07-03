function [L, l0, W] = regress_with_delay(u,x,d)
% [L, l0, W] = regress_with_delay(u,x,d)
%
% fits {L,l0,W} according to the following regression equation:
% u(t) = L*x(t-d) + l0 + w_t,  w_t ~ N(0,W)
%
% INPUTS: 
%   u and x are cell arrays, where u{i} and x{i} have the same number of
%   columns
%   d is an integer timestep delay

if isa(u,'double')
    u = {u};
    x = {x};
end

num_trials = numel(u);

sizeX1 = cellfun(@(x)(size(x,1)),x);
sizeX2 = cellfun(@(x)(size(x,2)),x);
sizeU1 = cellfun(@(u)(size(u,1)),u);
sizeU2 = cellfun(@(u)(size(u,2)),u);

% Check for appropriate input dimensionalities
if all(sizeX1==sizeU1) && numel(unique(sizeX2))==1  && numel(unique(sizeU2))==1 && ~all(sizeX2==sizeU2)
    % x and u have the same number of rows
    % all x cells have the same dimensionality
    % all u cells have the same dimensionality
    u = cellfun(@transpose,u,'uniformoutput',false);
    x = cellfun(@transpose,x,'uniformoutput',false);
elseif all(sizeX2==sizeU2) && numel(unique(sizeX1))==1  && numel(unique(sizeU1))==1
    % x and u have the same number of columns
    % all x cells have the same dimensionality
    % all u cells have the same dimensionality
    % do nothing
elseif all(sizeX1==sizeU1) && all(sizeX2==sizeU2)
    % Ambiguous which dimension is what in cell arrays
    % Proceed as if observations are rows of cell array matrices
else
    % error('invalid input cell array(s)');
end

[U, X] = cell2mat_forLaggedRegression(u, x, d);

% Append ones for bias term (i.e. to account for mean spike counts)
if isempty(X)
    % If no predictors/regressors/features are given, then this regression
    % will only return l0, i.e. the data mean: l0 = mean(U,2)
    X = ones(1,size(U,2));
    % --> X*X' = N
    % --> U*X' = sum(U,2)
    % --> l0 = sum(U,2)/N = mean(U,2)
else
    X(end+1,:) = 1;
end

Ll0 = (U*X')/(X*X');

l0 = Ll0(:,end);
L = Ll0(:,1:end-1);

residual = (U - Ll0*X);

N = size(residual,2);
W = (residual * residual')/N;

% Notice that l0 is close to but not exactly the mean firing rate of the
% data used (unless input predictors, X, are empty, in which case l0=mean(U,2))
% [l0 mean(U,2)]

end
