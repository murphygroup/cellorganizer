function [direct_neighbors, face_neighbors, face_ids] = compute_neigbhors(faces, mode)
% compute ordered neighbors for faces. 
% two mode: mode 1, just calulate direct_neighbors, without ordering them, default mode.
%           mode 2, compute ordered neighbors as well as diagonal
%           neighbors.

if nargin < 2
    mode = 1;
end

unique_inds = unique(faces(:));
nvert = numel(unique_inds);
nface = size(faces, 1);

circ_ind = @(x, n) rem(x + n - 1, n) + 1;
face_neighbors = cell(nvert, 1);
direct_neighbors = cell(nvert, 1);
face_ids = cell(nvert, 1);

for i = 1 : nface
    cur_face = faces(i, :);
    for j = 1 : numel(cur_face)
        face_ids{cur_face(j)}{end + 1} = i;
        direct_neighbors{cur_face(j)} = [direct_neighbors{cur_face(j)}, cur_face(circ_ind(j+1, 4)), cur_face(circ_ind(j-1, 4))];
    end
end

face_ids = cellfun(@(x) unique(cat(1, x{:})), face_ids, 'uniformoutput', false);
direct_neighbors = cellfun(@(x) unique(x), direct_neighbors, 'uniformoutput', false);
% if in mode 1, there is no need to compute ordered faces and neighbors,
% which takes time. 
if mode == 1
    return;
end    

% convert face neighbor to a ordered list
visit_mat = zeros(nvert, 1);
face_visit_mat = zeros(nface, 4);
neighbor_visit_mat = zeros(nvert, 1);
face_order_mat = zeros(nface, 4);

% associate edges with faces and also order faces if there is one with
% differet direction
% face_edge_mat = sparse(nvert, nvert);
face_edge_mat = spalloc(nvert, nvert, 20 * nvert);

for i = 1 : nface
    cur_face = faces(i, :);
    inds = [cur_face', circshift(cur_face', -1)];
    inds = sub2ind(size(face_edge_mat), inds(:, 1), inds(:,2 ));
    
    if any(face_edge_mat(inds) ~= 0)
        disp('The order is inverted!')
        cur_face = flip(cur_face);
        inds = [cur_face', circshift(cur_face', -1)];
        inds = sub2ind(size(face_edge_mat), inds(:, 1), inds(:,2 ));      
    end
    face_edge_mat(inds) = i;  
    face_order_mat(i, :) = cur_face;
end

cur_ind = 1;

% find ordered neighbors by travesal of the ordered faces. 
while ~all(visit_mat)
    cur_nb_list = [];
    cur_face_ids = face_ids{cur_ind};
    end_ind = [];
    while size(cur_face_ids > 0)
        if isempty(end_ind)
            if all(face_visit_mat == 0)
                % for the first visit
                cur_face_ind = 1;
            else
                % for a new visit
                visited_face_inds = find(face_visit_mat(cur_face_ids) == 1);
                cur_face_ind = visited_face_inds(1);
            end
            cur_face_id =  cur_face_ids(cur_face_ind);
            cur_verts = face_order_mat(cur_face_id, :);
            cur_loc = find(cur_verts == cur_ind);
            end_ind = cur_verts(circ_ind(cur_loc-1, 4));                
        else
            nb_face_inds = find(any(face_order_mat(cur_face_ids, :) == end_ind, 2));
            if isempty(nb_face_inds)
                break;
            end
            cur_face_ind = nb_face_inds(1);
            cur_face_id =  cur_face_ids(cur_face_ind);
        end
        
        cur_verts = face_order_mat(cur_face_id, :);
        cur_loc = find(cur_verts == cur_ind);
        cur_verts_shift = circshift(cur_verts, 1 - cur_loc);
        nbs = cur_verts_shift(2 : end);
        end_ind = nbs(end);         
        cur_nb_list = [cur_nb_list, nbs];
        face_visit_mat(cur_face_id) = 1;
        neighbor_visit_mat(nbs) = 1;
        cur_face_ids(cur_face_ind) = [];
    end
    cur_nb_list = unique(cur_nb_list, 'stable');    
    face_neighbors{cur_ind} = cur_nb_list;
    visit_mat(cur_ind) = 1;

    % move to another vertice 
    novisit_inds = find(visit_mat(cur_nb_list) == 0);
    if numel(novisit_inds) > 0
        cur_ind = cur_nb_list(novisit_inds(1));
    else
        novisit_inds = find(neighbor_visit_mat .* (1 - visit_mat) == 1);
        if numel(novisit_inds) > 0
            cur_ind = novisit_inds(1);
        end
    end
end        

% use ordered version of direct neibhors
direct_neighbors = cellfun(@(x) x(1 : 2 : end)', face_neighbors, 'uniformoutput', false);

end
