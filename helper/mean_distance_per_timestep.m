function d = mean_distance_per_timestep(X)

POS_IDX = 1:2;

distances = cellfun(@(x)(sqrt(sum(diff(x(POS_IDX,:),1,2).^2,1))),X,'uniformoutput',false);
d = mean(cell2mat(distances));