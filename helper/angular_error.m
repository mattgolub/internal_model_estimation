function degs = angular_error(v1,v2)
% Vectors are rows of v1 and v2.  Both inputs must be the same size
% Fixed bug that makes 180 degree errors become 0: 6/27/2014

if ~isequal(size(v1),size(v2))
    error('Inputs are not the same size')
end

error_is_counterclockwise = sign(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));

v1norm = sqrt(sum(v1.^2,2));
v2norm = sqrt(sum(v2.^2,2));

dot_prods = sum(v1.*v2,2);
degs = real(rad2deg(acos(dot_prods./(v1norm.*v2norm))));
degs(error_is_counterclockwise==1) = -degs(error_is_counterclockwise==1);