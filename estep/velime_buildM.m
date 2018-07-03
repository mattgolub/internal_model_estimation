function [M, M0] = velime_buildM(estParams)
% [M, M0] = velime_buildM(estParams)
% 
% Constructs concatenated linear system so the unrolling of all prior
% expectations can be computed with the single evaluation: X = M*C + M0
% C has columns c_t = [p_{t-tau}; v_{t-tau-1}; u_{t-tau+1}; ... ; u_{t+1}]
% and is computed in fastfmc_assemble_data.
%
% @ Matt Golub, 2014

A = estParams.A;
B = estParams.B;
b0 = estParams.b0;
TAU = estParams.TAU;
dt = estParams.dt;

vDim = size(A,1);
pDim = vDim;
xDim = pDim + vDim;
uDim = size(B,2);

% Build M0 = [0; b0; b0*dt; A*b0+b0; (A*b0+2*b0)*dt; ...]
M0 = zeros((TAU+1)*xDim,1);
M0(pDim+1:xDim,1) = b0;

M = zeros((TAU+1)*xDim,xDim + (TAU+1)*uDim);
M(1:pDim,1:pDim) = eye(pDim);
M(pDim+1:pDim+vDim,pDim+1:pDim+vDim+uDim) = [A B];
col = xDim + uDim;
for row = xDim:xDim:TAU*xDim
    prev_p_idx = row-xDim+1:row-vDim;
    prev_v_idx = row-vDim+1:row;
    
    cur_p_idx = row+1:row+pDim;
    cur_v_idx = row+pDim+1:row+xDim;
    cur_x_idx = row+1:row+xDim;
    
    % M
    prev_p_term = M(prev_p_idx,:);
    prev_v_term = M(prev_v_idx,:);
    M(cur_p_idx,:) = prev_p_term + prev_v_term*dt;
    M(cur_v_idx,:) = A*prev_v_term;
    M(cur_v_idx,col+1:col+uDim) = B;
    col = col + uDim;
    
    % M0
    prev_p_term = M0(prev_p_idx,1);
    prev_v_term = M0(prev_v_idx,1);
    M0(cur_x_idx,1) = [prev_p_term + prev_v_term*dt; A*prev_v_term + b0];
end

end