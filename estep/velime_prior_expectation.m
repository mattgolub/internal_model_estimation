function [E_X, COV_X, E_TARG] = velime_prior_expectation(C, estParams)
% [E_X, COV_X, E_TARG] = velime_prior_expectation(C, estParams)
% 
% Computes 
%   1) the prior distributions [x_{t-tau}^t,...,x_t^t] given
%   p_{t-tau-1}, p_{t-tau}, u_{t-tau+1}, ..., u_{t+1}
% 
%   2) the prior expectation of the target position, 
%   E(G_t | p_{t-tau-1}, p_{t-tau}, u_{t-tau+1}, ..., u_{t+1}).
% 
% @ Matt Golub, 2014.

COMPUTE_COV_X = nargout>=2;
COMPUTE_E_TARG = nargout==3;

TAU = estParams.TAU;

[M, M0] = velime_buildM(estParams);
E_X = bsxfun(@plus,M*C, M0);

% COV_X = cov([x_{t-tau}^t;...;x_{t+1}^t] | y_{t-tau}, u_{t-tau},...u_{t})
if COMPUTE_COV_X
    A = estParams.A;
    W_p = estParams.W_p;
    W_v = estParams.W_v;
    dt = estParams.dt;
    
    vDim = size(A,1);
    pDim = vDim;
    xDim = vDim + pDim;
    
    M1 = [eye(pDim) eye(vDim)*dt;
        zeros(pDim) A];
    W = blkdiag(diag(W_p),diag(W_v));
    
    COV_X = zeros((TAU+1)*xDim);
    
    [K, x_idx] = velime_x_index(TAU,xDim);
    
    % First compute block diagonal components:
    COV_X(1:pDim,1:pDim) = 0; % position feedback is observed (this term does NOT depend on W_p=0)
    COV_X(pDim+1:xDim,pDim+1:xDim) = diag(W_v); % first velocity in state is NOT observed
    for k = -TAU+1:0        
        xkm1_idx = x_idx(:,K==k-1);
        xk_idx = x_idx(:,K==k);
        COV_X(xk_idx,xk_idx) = M1*COV_X(xkm1_idx,xkm1_idx)*M1' + W;
    end
    
    % Now compute off diagonals
    for k1 = -TAU:0 % first velocity in state is NOT observed
        xk1_idx = x_idx(:,K==k1);
        for k2 = k1+1:0
            xk2_idx = x_idx(:,K==k2);
            xk2m1_idx = x_idx(:,K==k2-1);
            COV_X(xk1_idx,xk2_idx) = COV_X(xk1_idx,xk2m1_idx)*M1';
            COV_X(xk2_idx,xk1_idx) = COV_X(xk1_idx,xk2_idx)';
        end
    end
end

if COMPUTE_E_TARG
    % This computes E(G_t | y_{t-tau}^t, u_{t-tau},...u_{t})
    %
    % Note: this only makes sense if estParams was fit on EXACTLY the data
    % which went into C.
    
    xPosIdx = 1:pDim; % 1:2
    xVelIdx = pDim+1:xDim; % 3:4
    E_TARG = E_X(x_idx(xPosIdx,K==0),:) + bsxfun(@times,E_X(x_idx(xVelIdx,K==0),:),estParams.alpha);
end

end