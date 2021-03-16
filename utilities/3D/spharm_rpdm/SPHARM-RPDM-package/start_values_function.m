function [cartesian] = start_values_function(vertices, faces, neighbors)
% spherical parameterization initialation using the method from Brechbuhler
% & Gerig et al. 1995. 
    
pole_inds = [1, size(vertices, 1)];

% compute latitude
A1 = setup_matrix_A(vertices, neighbors, pole_inds);
b1 = setup_vector_b(vertices, neighbors, pole_inds);

latitudes = A1 \ b1;
latitudes = [0; latitudes; pi];

% compute longitude
A2 = modify_matrix_A(A1, neighbors, pole_inds);
b2 = setup_logitude_vector_b(vertices, neighbors, latitudes, pole_inds);

longitudes = A2 \ b2;
longitudes = [0; longitudes; 0];

cartesian = [sin(latitudes) .* cos(longitudes), sin(latitudes) .* sin(longitudes), cos(latitudes)];

end


function [A] = setup_matrix_A(vertices, neighbors, pole_inds)
% we assume first point north pole, last point sourth pole

nvert = size(vertices, 1);
A = sparse(nvert, nvert);

inds = arrayfun(@(ind) (neighbors{ind} -1) .* nvert  + ind, 1 : nvert, 'uniformoutput', false);
inds = cat(1, inds{:});
A(inds) = -1;

n_neighbor = cellfun(@numel, neighbors);

% unique_pole_neighbors = unique([neighbors{pole_inds}]);
A(:, pole_inds) = 0;

A(speye(nvert) == 1) = n_neighbor;
A(pole_inds, :) = [];
A(:, pole_inds) = [];

end


function [b] = setup_vector_b(vertices, neighbors, pole_inds)

nvert = size(vertices, 1);
n = nvert - 2;
b = zeros(n, 1);
b(neighbors{pole_inds(2)} - 1) = pi;

end


function [A] = modify_matrix_A(A, neighbors, pole_inds)

unique_pole_neighbors = unique([neighbors{pole_inds}]);
pole_nb_inds = sub2ind(size(A), unique_pole_neighbors-1, unique_pole_neighbors-1);
A(pole_nb_inds) = A(pole_nb_inds) - 1;
A(1, 1) = A(1, 1) + 2;

end


function [b] = setup_logitude_vector_b(vertices, neighbors, theta_lati, pole_inds)

nvert = size(vertices, 1);
n = nvert - 2;

b = zeros(n, 1);

% north pole
previous = 1;
here = neighbors{previous}(1);
maximum = 0;
max_ind = 1;

while here ~= pole_inds(2)
    cur_neighbors = neighbors{here};
    maximum = 0;
    max_ind = 1;
    for i = 1 : numel(cur_neighbors)
        neighbor = cur_neighbors(i);
        if theta_lati(neighbor) > maximum
            maximum = theta_lati(neighbor);
            max_ind = neighbor;
            nextpos = vertices(neighbor, :);
            next = neighbor;
        end
        if neighbor == previous
            prevpos = vertices(neighbor, :);
            previous = neighbor;
        end 
    end  
    previous_ind = find(cur_neighbors == previous);
    cur_neighbors_shift = circshift(cur_neighbors, 1 - previous_ind);
    for i = 2 : numel(cur_neighbors_shift)
        if cur_neighbors_shift(i) == next
            break;
        end
        neighbor = cur_neighbors_shift(i);
        b(neighbor - 1) = b(neighbor - 1) + 2 * pi;
        b(here - 1) = b(here - 1) - 2 * pi;
    end
    previous = here;
    here = next;
end

end

