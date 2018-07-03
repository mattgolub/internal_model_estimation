function [K, x_idx, var_names] = velime_x_index(TAU,xDim)
% Example for TAU=3, xDim=4:
% K = [-3    -2    -1     0]
% x_idx =
%      1     5     9    13
%      2     6    10    14
%      3     7    11    15
%      4     8    12    16

K = -TAU:0;
x_idx = reshape(1:(TAU+1)*xDim,xDim,TAU+1);

if nargout==3
    var_names = cell(2,TAU+1);
    k = min(K);
    for i=1:TAU
        str = num2str(k);
        var_names{1,i} = ['p_{t' str '}'];
        var_names{2,i} = ['v_{t' str '}'];
        k = k+1;
    end
    var_names{1,i+1} = 'p_t';
    var_names{2,i+1} = 'v_t';
end