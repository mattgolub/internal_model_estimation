function estParams = velime_mstep_alphaParams(G, E_X_posterior, COV_X_posterior, const)
% estParams = velime_mstep_alphaParams(G, E_X_posterior, COV_X_posterior, const)
%
% @ Matt Golub, 2014.

pt_idx = const.x_pt_idx;
vt_idx = const.x_vt_idx;
xt_idx = const.x_idx(:,end);
gDim = const.gDim;

E_xt = E_X_posterior(xt_idx,:);
E_pt = E_X_posterior(pt_idx,:);
E_vt = E_X_posterior(vt_idx,:);

COV_xt = COV_X_posterior(xt_idx,xt_idx,:);
COV_ptpt = COV_X_posterior(pt_idx,pt_idx,:);
COV_ptvt = COV_X_posterior(pt_idx,vt_idx,:);
COV_vtvt = COV_X_posterior(vt_idx,vt_idx,:);

T = const.T;
% Indices of the diagonal elements of a xDim x xDim subblock of COV_X_posterior
SLICE_DIAG_IDX = 1:(gDim+1):(gDim^2);

% Diagonal indices of each slice of a 3-D vector for trace covaraiance computations
DIAG_IDX = bsxfun(@plus,SLICE_DIAG_IDX',0:gDim^2:T*((gDim^2))-1);

trace_E_pp = sum(COV_ptpt(DIAG_IDX) + E_pt.*E_pt);
trace_E_pv = sum(COV_ptvt(DIAG_IDX) + E_pt.*E_vt);
trace_E_vv = sum(COV_vtvt(DIAG_IDX) + E_vt.*E_vt);
trace_E_Gv = sum(E_vt.*G);

alpha = (trace_E_Gv - trace_E_pv)./trace_E_vv;

% Constrain alphas to be non-negative
alpha = max(alpha,0);

term1 = sum(G.*G);
term2 = -2*sum(G.*(E_pt + bsxfun(@times,E_vt,alpha)));
term3 = trace_E_pp + 2*alpha.*trace_E_pv + (alpha.^2).*trace_E_vv;
r = 1/(T*gDim)*sum(term1 + term2 + term3);
R = r*ones(gDim,1);

%% Slow, more readable version of the fast code above:
% T = size(E_X_posterior,2);
% alpha = zeros(1,T);
% r = nan(1,T);
% for t = 1:T
%     trace_E_pv = trace(COV_ptvt(:,:,t)) + E_pt(:,t)'*E_vt(:,t);
%     trace_E_vv = trace(COV_vtvt(:,:,t)) + E_vt(:,t)'*E_vt(:,t);
%     trace_E_Gv = E_vt(:,t)'*G(:,t);
%     
%     alpha(t) = (trace_E_Gv - trace_E_pv)/trace_E_vv;
%     alpha(t) = max(alpha(t),0);
%     C_t = [eye(2) eye(2)*alpha(t)];
%     
%     r(t) = G(:,t)'*G(:,t) - 2*G(:,t)'*C_t*E_xt(:,t) + trace(C_t*(COV_xt(:,:,t) + E_xt(:,t)*E_xt(:,t)')*C_t');
% end
% r = 1/(T*gDim) * sum(r);
% R = r*ones(gDim,1);

estParams.alpha = alpha;
estParams.R = R;