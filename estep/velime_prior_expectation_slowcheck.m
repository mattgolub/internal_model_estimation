function E_X = velime_prior_expectation_slowcheck(C, estParams)
% E_X = velime_prior_expectation_slowcheck(C, estParams)
%
% Computes E[x_{t-tau}^t,...,x_t^t | p_{t-tau-1}, p_{t-tau}, u_{t-tau+1}, ..., u_{t+1}]
% This computes recursively as a means to sanity check the often faster
% solution in velime_prior_expectation.m
%
% @ Matt Golub, 2014.

A = estParams.A;
B = estParams.B;
b0 = estParams.b0;
TAU = estParams.TAU;
dt = estParams.dt;

vDim = size(A,1);
pDim = vDim;
xDim = pDim + vDim;
EXDim = (TAU+1)*xDim;
uDim = size(B,2);

prev_p = C(1:2,:);
prev_v = C(3:4,:);

p_idx = 1:pDim; % 1:2
v_idx = (1+pDim):(1+xDim-1); % 3:4
u_idx = (xDim+1):(xDim + uDim);

E_X(p_idx,:) = prev_p;
E_X(v_idx,:) = bsxfun(@plus,A*prev_v + B*C(u_idx,:), b0);

u_idx = u_idx + uDim;
for x_idx = (xDim+1):xDim:EXDim
    p_idx = x_idx:(x_idx+pDim-1);
    v_idx = (x_idx + pDim):(x_idx+xDim-1);
    
    prev_p = E_X(p_idx-xDim,:);
    prev_v = E_X(v_idx-xDim,:);
    
    E_X(p_idx,:) = prev_p + prev_v*dt;
    E_X(v_idx,:) = bsxfun(@plus,A*prev_v + B*C(u_idx,:), b0);
    
    u_idx = u_idx + uDim;
end