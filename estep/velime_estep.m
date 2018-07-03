function [LLi, E_X_posterior, COV_X_posterior]  = velime_estep(C, G, estParams, const)
% [LLi, E_X_posterior, COV_X_posterior]  = velime_estep(C, G, estParams, const)
%
% Computes the posterior distributions of [x_{t-tau}^t,...,x_t^t] given
% p_{t-tau-1}, p_{t-tau}, u_{t-tau+1}, ..., u_{t+1}, G_t.
%
% @ Matt Golub, 2014.

R = diag(estParams.R);
alpha = estParams.alpha;

gDim = const.gDim;
llconst = gDim*log(2*pi);
half = 1/2;

x_pt_idx = const.x_pt_idx;
x_vt_idx = const.x_vt_idx;
xt_idx = [x_pt_idx x_vt_idx];

% Step 1: Put together E(x_{t-tau}^t,...,x_t^t, u_t | x_{t-tau}^t, u_{t-tau+1},...u_{t-1})
% and cov(x_{t-tau}^t,...,x_t^t, u_t | x_{t-tau}^t, u_{t-tau+1},...u_{t-1})
[E_X, SIGMA11, E_TARG] = velime_prior_expectation(C, estParams);
EXDim = size(E_X,1);

x2_minus_mu2 = G - E_TARG;

% Sigma_22 = cov(G_t | x_-{tau}^t, u_{t-tau+1},...u_t)
% SIGMA12 = cov([x_-{tau+1}^t;...;x_0^t], G_t | x_-{tau}^t, u_{t-tau+1},...u_{t-1})

T = size(E_X,2);
I = eye(const.gDim);

if const.gDim==2 
    % Fast code.
    % Leverages efficient computation of determinants of 2x2 matrices.
    
    SIGMA_pp = SIGMA11(x_pt_idx,x_pt_idx); % cov(p_t^t)
    SIGMA_pv = SIGMA11(x_pt_idx,x_vt_idx); % cov(p_t^t,v_t^t)
    SIGMA_vp = SIGMA11(x_vt_idx,x_pt_idx); % cov(v_t^t,p_t^t)
    SIGMA_vv = SIGMA11(x_vt_idx,x_vt_idx); % cov(v_t^t)
    SIGMA_Xp = SIGMA11(:,x_pt_idx); % cov(X,p_t^t)
    SIGMA_Xv = SIGMA11(:,x_vt_idx); % cov(X,v_t^t)
    
    % Build up terms to SIGMA22
    term1 = repmat(SIGMA_pp(:),1,T);
    term2 = bsxfun(@times,repmat(SIGMA_pv(:)+SIGMA_vp(:),1,T),alpha);
    term3 = bsxfun(@times,repmat(SIGMA_vv(:),1,T),alpha.^2);
    
    SIGMA22s = bsxfun(@plus,term1 + term2 + term3, R(:));
    % To recover SIGMA22: reshape(SIGMA22s(:,1),gDim,gDim)
    
    detSIGMA22s = SIGMA22s(1,:).*SIGMA22s(4,:)-SIGMA22s(2,:).*SIGMA22s(3,:);
    logdetSIGMA22s = log(detSIGMA22s);
    
    inv_term = [SIGMA22s(4,:); -SIGMA22s(2,:); -SIGMA22s(3,:); SIGMA22s(1,:)];
    invSIGMA22s = bsxfun(@times,inv_term,1./detSIGMA22s);
    
    XcXc = [x2_minus_mu2(1,:).^2
        x2_minus_mu2(1,:).*x2_minus_mu2(2,:)
        x2_minus_mu2(1,:).*x2_minus_mu2(2,:)
        x2_minus_mu2(2,:).^2];
    LLi = -half*(T*llconst + sum(logdetSIGMA22s + sum(XcXc.*invSIGMA22s,1)));
    
    % Each page, e.g., invSIGMA22_paged(:,:,t), is an invSIGMA22, 
    invSIGMA22_page = reshape(invSIGMA22s,gDim,gDim,T);
    
    SIGMA12s = repmat(SIGMA_Xp(:),1,T) + bsxfun(@times,repmat(SIGMA_Xv(:),1,T),alpha);
    
    SIGMA12_paged = reshape(SIGMA12s,EXDim,gDim,T);
    SIGMA12_iSIGMA22_paged = multiprod(SIGMA12_paged,invSIGMA22_page);
    SIGMA12_iSIGMA22_SIGMA21_paged = multiprod(SIGMA12_iSIGMA22_paged,multitransp(SIGMA12_paged));
    COV_X_posterior = baxfun(@minus,SIGMA11,SIGMA12_iSIGMA22_SIGMA21_paged);
    
    if nargout>1 % don't do this if LL is only desired output
        SIGMA12_iSIGMA22_X2_minus_mu2_paged = multiprod(SIGMA12_iSIGMA22_paged,reshape(x2_minus_mu2,[gDim 1 T]));
        E_X_posterior = E_X + squeeze(SIGMA12_iSIGMA22_X2_minus_mu2_paged);
    end
else
    % Slow, readable code. 
    % Works when determinants must be taken of square matrices larger than 
    % 2x2. Reproduces results from "fast code" above when gDim==2.
    
    LLi = 0;
    E_X_posterior = nan(EXDim,T);
    COV_X_posterior = nan(EXDim,EXDim,T);

    for t = 1:T
        C_t = [I alpha(t)*I];
        SIGMA22 = C_t*SIGMA11(xt_idx,xt_idx)*C_t' + R;
        [invSIGMA22, logdetSIGMA22] = invAndLogDet(SIGMA22);
        
        LLi = LLi - half*(llconst + logdetSIGMA22 + sum(sum((x2_minus_mu2(:,t)*x2_minus_mu2(:,t)').*invSIGMA22)));
        % LLi = LLi + logmvnpdf(G(:,t),E_TARG(:,t),SIGMA22);
        
        SIGMA12 = SIGMA11(:,xt_idx)*C_t';
        SIGMA12_iSIGMA22 = SIGMA12*invSIGMA22;
        COV_X_posterior(:,:,t) = SIGMA11 - SIGMA12_iSIGMA22*SIGMA12';
        
        % Find conditional expectation of latent variables given observed variables
        if nargout>1 % don't do this if LL is only desired output
            E_X_posterior(:,t) = E_X(:,t) + SIGMA12_iSIGMA22*(x2_minus_mu2(:,t));
        end
    end
end
end