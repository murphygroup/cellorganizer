function [A] = calc_jacob_matrix_func_1(vertices, faces, neighbors, face_ids, m_x, c_hat, active, param)
% calculate Jacobian of the system using finite differences

nface = size(faces, 1);
nvert = size(vertices, 1);
desired_area = 4 * pi / nface;

% inequality constraints
active = active(1 : end - 1);
active_faces = ceil(active / 4);
active_angle_inds = rem(active + 3, 4) + 1;

A = sparse(nface + numel(active),  nvert * 3);
c_hat = [c_hat(1 : nface-1); 0; c_hat(nface:end)];
m_x_try = m_x;

for col = 1 : nvert
    for i = 1 : 3
        m_x_try(col, i) = m_x(col, i) + param.delta;
        col_0 = col;
        col_elm_ind = (col - 1) * 3 + i;
        m_x_try(col_0, :) = m_x_try(col_0, :) ./ sqrt(sum(m_x_try(col_0, :) .^ 2));
        
        face_ids_i = face_ids{col};
        if ~false && any(ismember(6285, face_ids_i))
            flag = 1;
        end
        [areas_i, sines_i] = spher_area4_function(m_x_try, faces(face_ids_i, :));
        % area_c = areas_i - desired_area;
        A(face_ids_i, col_elm_ind) =  (areas_i - desired_area - c_hat(face_ids_i)) ./ param.delta;
        
        % see if any face is active
        [Ca, ia, ib] = intersect_function(face_ids_i, active_faces);
        active_inds = active(ib);
        active_f_inds = active_angle_inds(ib)';
        if ~isempty(Ca)
            sines_act_i = sines_i(sub2ind(size(sines_i), ia, active_f_inds)) - c_hat(ib + nface - 1);
            A(nface + ib, col_elm_ind) = sines_act_i ./ param.delta;
        end
        m_x_try(col_0, :) = m_x(col_0, :);
    end
end
A(nface, :) = [];

end


function [Ca, iaa, iba] = intersect_function(A, B)
% implement intersect function with keep of all same element in ia or ib
% A, B are row vectors
A = reshape(A, [], 1);
B = reshape(B, [], 1);

Lia = ismember(A, B);
Ca = A(Lia);
iaa = {};
iba = {};
for i = 1 : numel(Ca)
    iba{i} = find(B == Ca(i));
    iaa{i} = repmat(find(A == Ca(i)), numel(iba{i}), 1);
end

iaa = cat(1, iaa{:});
iba = cat(1, iba{:});

end