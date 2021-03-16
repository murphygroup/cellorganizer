function [ keepinds, distance_matrix_trimmed ] = get_complete_distance_matrix( distance_matrix )

distance_matrix(eye(size(distance_matrix))>0) = nan;

keepinds = find(~all(isnan(distance_matrix),2));

distance_matrix = distance_matrix(keepinds,keepinds);
distance_matrix(eye(size(distance_matrix))>0) = 0;

nancounts = sum(isnan(distance_matrix),2);

while any(nancounts)
    [~, ind] = max(nancounts);

    distance_matrix(ind,:) = [];
    distance_matrix(:, ind) = [];
    keepinds(ind) = [];

    nancounts = sum(isnan(distance_matrix),2);
end

distance_matrix_trimmed = distance_matrix;

end

