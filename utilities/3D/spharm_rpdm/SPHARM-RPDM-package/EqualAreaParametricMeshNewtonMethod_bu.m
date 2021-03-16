function [m_x, cost_mat, is_success] = EqualAreaParametricMeshNewtonMethod(vertices, faces, m_x, par)
% Input: vertices, nvert * 3 
%        faces: quad faces, nface * 4
%        bim: binary volume of the image, if not provided, use vertices and
%        faces to generate it. 
%        par: parameter for the optimization
% Output: m_x: spherical parameterization
%         cost_mat: lagrangian cost for the optimization. 
% 
% perform spherical parameterization optimization, the method is basically
% modified from SPHARM-PDM EqualAreaParametricMeshNewtonIterator.cxx
% which is based on Brechbuhler & Gerig et al. 1995. 
% 
% Improvement: 1. uses analytic Jacobian for constraints (both area and
% sines), which are much accurate and faster for computing. 
% 2. use the initialization frame from SPHARM-MAT, which is based on the
% papaer Shen & Makedon 2006. Using this initialization framework, the ,
% there is less distortion in faces, the program can converge faster, or handle
% difficult cases that the old frame cannot. 
%
% Author: Xiongtao Ruan (xruan@cs.cmu.edu)
% Date: March, 2018
% 
% 03/13/2018 add return of objective function, and use it to decide whether
% save best m_x, rather than use cost (which is not value of objective function).
% 03/25/2018 change the returnof objective function to largarigian, which
% can better reflect the quality of parameterization. 
% 04/06/2018 add another criterion for convergence, that is, the largr
% difference of two adjacency iteration smaller than a threshold (1e-7).
% Also set a threshold for the update of largr_best, to make sure it's not
% random fluctuation. 
% 04/09/2018 first use lsqminnorm and see if the solution is constrained,
% if not use box constrained admm

debug = false;

% Initialize parameter
default_par.max_active = 500;
default_par.print_itn  = -2;
% default_par.delta      = 3e-7;
default_par.delta      = 1e-8;
default_par.constr_tol = 1e-3;
default_par.line_tol   = 1e-5;
default_par.ineq_low   = 1e-7;
default_par.ineq_init  = 1e-2;
default_par.ineq_final = 1e-6;
default_par.ineq_slack = 2.0;
default_par.newton_tol = 1e-4;
default_par.rho_init   = 1;
default_par.c0rho      = 1;
default_par.c1rho      = 0.25;
default_par.c2rho      = 0.012;
default_par.rho_limit  = 3e-2;
default_par.step_small = 0.5;
default_par.step_large = 1.0;
default_par.cost_tol = 1e-7;  % the minimum cost if smaller than that cost, the optimization finished. 
default_par.largr_tol = 1e-7;  % if the absolute difference of largragian between two iterations samller than threshold, the optimization finished. 
default_par.max_iter = 1500; % max iteration 
default_par.jacobi_cond_tol = 1e10; % max conditional number 
default_par.initialization_method = 'default';  
default_par.lsq_bound = 100; % bound of solutions of least squares in the update.  

if nargin < 4
    par = struct();
end
par = process_options_structure(default_par, par);

jacobi_cond_tol = par.jacobi_cond_tol;
lsq_bound = par.lsq_bound;

J_sines_func = cell(4, 1);
for i = 1 : 4
    function_name = sprintf('calc_sine_%d_jacobian_sphere', i);
    J_sines_func{i} = str2func(function_name);
end

% initialization
m_rho = par.rho_init;
m_alpha_step = 1.0;

no_activation = -1;
nvert = size(vertices, 1);
nface = size(faces, 1);
% compute ordered neighbors of all vertices. 
[neighbors, ~, face_ids] = compute_neigbhors(faces, 1);
total_nb_num = sum(cellfun(@numel, neighbors));

% initialization method
if nargin < 3 || isempty(m_x)
    m_x = setup_init_parameterization(vertices, faces, par);
end

gCG = m_x;
m_x_try = m_x;
hCG = zeros(nvert, 3);
m_dx = hCG;

[equal_mat, constr_ineq] = constraints(m_x, faces);
% set last face unconstrained
c_hat = equal_mat(1:end-1);
ineq_high = par.ineq_init;
inequal_mat = reshape(constr_ineq', [], 1);
[c_hat, active, activity, n_active] = activate_initial(c_hat, inequal_mat, ineq_high, no_activation);

MaxIter = par.max_iter;
% may change MaxIter if nvert is very large
MaxIter = max(MaxIter, ceil(nvert / 50));
best_largr = 100000;
largr_last = best_largr;
best_largr_bu = best_largr;
best_largr_diff = 10;
m_x_best = m_x;
cost_mat = zeros(MaxIter, 3);
obj_val = 1000000;
lambda_t_c_hat = 1e5;
% check whether the jacobi matrix ill 
is_ill_ever = true;
last_ill_iter = 1;
ill_dist = 5;

% set up optimization options for lsqlin
optm_options = optimoptions('lsqlin','Algorithm', 'trust-region-reflective', 'Display','none', 'FunctionTolerance', 1e-8);


% iterate to update
for iter = 1 : MaxIter
    tic;
    fprintf(" %3d", iter);   
    
    % jacobi_mat_1 = calc_jacob_matrix_func_1(vertices, faces, neighbors, face_ids, m_x, c_hat, active, par);
    [jacobi_mat, is_jacobi_ill] = calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, c_hat, J_sines_func, par);
    % I_sparse = speye(size(jacobi_mat, 1));
    if is_jacobi_ill || is_ill_ever || abs(lambda_t_c_hat) > 100
        if sprank(jacobi_mat) < size(jacobi_mat, 1) || condest(jacobi_mat * jacobi_mat') > jacobi_cond_tol
            jacobi_ill = true;
            last_ill_iter = iter;
        else
            jacobi_ill = false;
        end
        if iter - last_ill_iter > ill_dist
            is_ill_ever = false;
        end 
    else
        jacobi_ill = false;
    end
    if jacobi_ill
        % m_newton_dir = jacobi_mat \ c_hat;
        % later will use lsqminnorm
        % m_newton_dir = lsqr(jacobi_mat,  c_hat, 1e-6, 2000);
        m_newton_dir = lsqminnorm(jacobi_mat, c_hat);
        if any(abs(m_newton_dir) > lsq_bound)
            m_newton_dir = box_constrained_lsq_admm(jacobi_mat, c_hat, lsq_bound, -lsq_bound);
        end
    else
        m_newton_dir = jacobi_mat' * ((jacobi_mat * jacobi_mat') \ c_hat);
    end
    
    m_dx = reshape(m_newton_dir, 3, [])';
    
    act = no_activation; act_keep = no_activation;
    newtonStep = -1;
    for i = 1 : 22
        [c_hat_l, constr_ineq, badness, min_ineq, m_x_try, act] = step_and_check(faces, newtonStep, m_x, m_dx, active, n_active, c_hat, ineq_high);
        if act == no_activation && badness <= par.newton_tol
            break;
        end
        newtonStep = newtonStep / 2;
        act_keep = act;
    end
    
    fprintf(" %6.0e %6d", newtonStep, act_keep);
    m_x = spher_step(newtonStep, m_x, m_dx);
    c_hat_temp = c_hat;
    c_hat = c_hat_l;
    c_hat_l = c_hat_temp;
    c_hat_r = c_hat;

    grad = calc_gradient_func(m_x, neighbors);
    grad = reshape(grad', [], 1);

    [jacobi_mat, is_jacobi_ill] = calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, c_hat, J_sines_func, par);
    % jacobi_mat_1 = calc_jacob_matrix_func_1(vertices, faces, neighbors, face_ids, m_x, c_hat, active, par);    
    if is_jacobi_ill || is_ill_ever || abs(lambda_t_c_hat) > 100
        if sprank(jacobi_mat) < size(jacobi_mat, 1) || condest(jacobi_mat * jacobi_mat') > jacobi_cond_tol
            jacobi_ill = true;
            last_ill_iter = iter;
        else
            jacobi_ill = false;
        end
        if iter - last_ill_iter > ill_dist
            is_ill_ever = false;
        end 
    else
        jacobi_ill = false;
    end

    if jacobi_ill
        % m_lambda = lsqr(jacobi_mat',  grad, 1e-6, 2000);
        m_lambda =  lsqminnorm(jacobi_mat', grad);
        if any(abs(m_lambda) > lsq_bound)
            m_lambda = box_constrained_lsq_admm(jacobi_mat', grad, lsq_bound, -lsq_bound);
        end
		% ub = 100 * ones(size(jacobi_mat, 1), 1);
		% lb = -100 * ones(size(jacobi_mat, 1), 1);
		% m_lambda = lsqlin(jacobi_mat', grad, [], [], [], [], lb, ub, [], optm_options);
    else
        m_lambda = (jacobi_mat * jacobi_mat') \ (jacobi_mat * grad);
    end        
    if max(abs(m_lambda)) > 1e5
        is_ill_ever = true;
        last_ill_iter = iter;
    end
    inactivated = 0;
    for i = n_active : -1 : 1
        if m_lambda(nface - 1 + i) < 0 && c_hat(nface - 1 + i) > -par.ineq_low
            [c_hat, flag, activity, active, n_active] = inactivate(nface, c_hat, activity, active, n_active, i);
            inactivated = inactivated + flag;
        end
    end

    if inactivated > 0
        [c_hat, active, activity, n_active, flag] = activate(act_keep, m_x_try, faces, c_hat, activity, active, n_active, "in_out");
        toc;
        continue;   
    end

    gradY = jacobi_mat' * m_lambda;
    % numerator = 0;
    % denominator = 0;
    % norm2gradZ = 0;
    gCG = reshape(gCG', [], 1);

    % conjugate gradient
    gradZ = gradY - grad;
    numerator = (gradZ -gCG)' * gradZ;
    denominator = gCG' * gCG;
    norm2gradZ = gradZ' * gradZ;
    gCG = gradZ;
    gamma = numerator / denominator;
    hCG = reshape(gCG, 3, [])' + gamma * hCG;
    m_dx = hCG;
    fprintf(" %6.3f %8.1e", gamma, sqrt(norm2gradZ));

    [m_alpha_step, c_sqr_sum, activated, c_hat, c_hat_l, c_hat_r, active, activity, n_active, min_ineq, m_x, m_x_try, m_rho, largr, obj_val, lambda_t_c_hat] = line_search(m_alpha_step, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, activity, n_active, c_hat, c_hat_l, c_hat_r, ineq_high, obj_val, par);
    [c_hat, active, activity, n_active, flag] = activate(act_keep, m_x_try, faces, c_hat, activity, active, n_active, "newton");
    % activated = activated + flag;
    
    % if activated == 0
    %     last_complete = iter;
    % end
    if -min_ineq * par.ineq_slack < ineq_high
        ineq_high = -min_ineq * par.ineq_slack;
    end
    if par.ineq_final > ineq_high
        ineq_high = par.ineq_final;
    end
    cost = sqrt(c_sqr_sum) + norm2gradZ;
    cost_mat(iter, :) = [largr, obj_val, cost];
    largr_diff = largr - largr_last;
    largr_last = largr;
    fprintf(" %10.3e %10.3e %10.3e %10.3e\n", sqrt(c_sqr_sum), cost, ineq_high, largr_diff);
    % only use the best largr and also not consider random flucturation.
    if largr > 0 && largr < best_largr        
        if abs(largr_diff) < 1e-1
            m_x_best = m_x;
            best_largr = largr;     
            best_largr_diff = largr_diff;
        elseif best_largr_bu > largr
            m_x_best_bu = m_x;
            best_largr_bu = largr;
        end
    end
    if cost < par.cost_tol || (largr > 0 && abs(largr_diff) < par.largr_tol)
        cost_mat(iter + 1 : end, :) = [];
        m_x_best = m_x;
        break;      
    end
    toc;

    if debug && rem(iter, 1000) == 0
        save('workspace.mat')
    end
end

if cost > 1e-2 || abs(lambda_t_c_hat) >  1e-2
    is_success = false;
else
    is_success = true;
end
if abs(best_largr_diff) > 1e-1
    is_success = false;
    m_x = m_x_best_bu;    
else
    m_x = m_x_best;
end

disp('Optimization complete!');

if debug
    save('workspace.mat')
end

end


function [obj] = goal_func(m_x, neighbors, total_nb_num)

nvert = size(m_x, 1);
nbsum_vec = zeros(nvert, 3);
for i = 1 : nvert
    nbsum_vec(i, :) = sum(m_x(neighbors{i}, :));    
end
% obj = sum(m_x .* nbsum_vec, 2);
% obj = 0.5 * sum(cellfun(@numel, neighbors) - obj(:));
obj = 0.5 * (total_nb_num - sum((m_x(:) .* nbsum_vec(:))));

end


function [grad] = calc_gradient_func(m_x, neighbors)
% calculate gradient of the objective function

nvert = size(m_x, 1);
nbsum_vec = zeros(nvert, 3);

for i = 1 : nvert
    nbsum_vec(i, :) = sum(m_x(neighbors{i}, :));    
end

grad = sum(m_x .* nbsum_vec , 2) .* m_x - nbsum_vec;

end


function [A, is_A_ill] = calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, c_hat, J_sines_func, param)
% calculate Jacobian of the system using analytic solution from symbolic
% computing, which is much faster than finite difference
% sometimes there are numerical error, which cause very large values or NaN
% for some elements, then we use finite difference to deal with these
% cases.
% 03/25/2018 fix bug for no active, c_hat for last unconstrained face. 

nface = size(faces, 1);
nvert = size(vertices, 1);
desired_area = 4 * pi / nface;
clip_val = 1e4;
is_A_ill = false;

if n_active > 0
    % inequality constraints
    active = active(1 : end - 1);
    active_faces = ceil(active / 4);
    active_angle_inds = rem(active + 3, 4) + 1;
    % A = zeros(nface + numel(active),  nvert * 3);
    circ_ind = @(x, n) rem(x + n - 1, n) + 1;
else
    % pad a value for the last unconstrained face
    c_hat = [c_hat; 0];
end

MaxIndex = 30 * nvert;
index_val_mat = zeros(MaxIndex, 3);
counter = 1;

% coordinate in a face, 4 * 3
n_num = 12;
m_x = m_x';

for i = 1 : nface
    cur_face = faces(i, :);
    col_elm_ind = (cur_face - 1) * 3 + (1:3)'; 
    col_elm_ind = col_elm_ind(:);
    cmx = m_x(:, cur_face);
    % J_area = J_area_func(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4));
    J_area = calc_area_jacobian_sphere(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4));
    % sometime when c_hat changes, the area may changes a lot, which makes
    % jacobian NaN, to avoid it, we use finite difference to compute the
    % jacobian for this face
    % we also deal with the same case when the jacobian is very large,
    % here we set a threshold of 1e4. 
    J_area = J_area(:);
    if any(isnan(J_area)) || any(abs(J_area) > clip_val)
        is_A_ill = true;
        for j = 1 : numel(col_elm_ind)
            % some are good
            if ~isnan(J_area(j)) && abs(J_area(j)) < clip_val
                continue;
            end
            cmxt = cmx;
            cmxt(j) = cmxt(j) + param.delta;
            cmxt = cmxt ./ sqrt(sum(cmxt .^ 2));
            [area_j, ~] = spher_area4_function(cmxt', [1, 2, 3, 4]);
            J_area_j = (area_j - desired_area - c_hat(i)) ./ param.delta;
            % clip the gradient to [-1e4, 1e4];
            J_area(j) = (J_area_j >= clip_val) .* clip_val + (J_area_j <= -clip_val) .* (-clip_val) + (abs(J_area_j) < clip_val) .* J_area_j;
        end
    end
    % A(i, col_elm_ind(:)) = J_area(:);
    counter_old = counter;
    counter = counter + n_num;
    index_val_mat(counter_old : counter - 1, :) = [ones(n_num, 1) .* i, col_elm_ind,  J_area];

    if n_active > 0 && any(active_faces == i)
        active_inds = find(active_faces == i);
        angle_inds = active_angle_inds(active_inds);
        for j = 1 : numel(angle_inds)
            active_ind = active_inds(j);
            angle_ind = angle_inds(j);
            cmxt = cmx;
            cmxt(:, circ_ind(angle_ind + 2, 4)) = [];            
            J_sine_func = J_sines_func{angle_ind};
            J_sines = J_sine_func(cmxt(1, 1), cmxt(1, 2), cmxt(1, 3), cmxt(2, 1), cmxt(2, 2), cmxt(2, 3), cmxt(3, 1), cmxt(3, 2), cmxt(3, 3));
            % A(nface + active_ind, col_elm_ind(:)) = J_sines(:);
            counter_old = counter;
            counter = counter + n_num;     
            index_val_mat(counter_old : counter - 1, :) = [ones(n_num, 1) .* (nface + active_ind), col_elm_ind,  J_sines(:)];
        end   
    end  
end

% put pairs values to the matrix
index_val_mat = index_val_mat(1 : counter-1, :);
A = sparse(index_val_mat(:, 1), index_val_mat(:, 2), index_val_mat(:, 3), nface + n_active,  nvert * 3);
A(nface, :) = [];

end


function [c_hat, active, activity, n_active, flag] = activate(act, m_x_try, faces, c_hat, activity, active, n_active, info)
no_activation = -1;
flag = 0;
if act == no_activation
    return;
end
fprintf("\n]]] %s : constraint %d -> level %d ", info, act, activity(act) + 1);
activity(act) = activity(act) + 1;
if activity(act) ~= 3
    return;
end
activity(act) = 5;

nface = size(faces, 1);
if n_active == 0
    i = 1;
else
    i = find(active > act, 1, 'first');
    if isempty(i)
        i = n_active + 1;
    end
end
active = [active(1 : i - 1), act, active(i:end)];
c_hat = [c_hat(1 : nface - 1 + i - 1); one_inequality(m_x_try, faces, act); c_hat(nface - 1 + i : end)];
n_active = n_active + 1;

fprintf(">>> %s activates constraint %d (pos %d); now active: %d\n                                                 ",  ...
    info, act, i, n_active);

flag = 1;
    
end


function [c_hat, flag, activity, active, n_active] = inactivate(nface, c_hat, activity, active, n_active, pos)

flag = 0;
fprintf('\n[[[ constraint %d -> level %d ;          ', active(pos), activity(active(pos)) - 1);

activity(active(pos)) = activity(active(pos)) - 1;
if activity(active(pos)) ~= 3
    return;
end

activity(active(pos)) = 2;
n_active = n_active - 1;
fprintf('<<< inactivate constraint %d at position %d; now active: %d\n', active(pos), pos, n_active);

c_hat(nface - 1 + pos) = [];
active(pos) = [];
flag = 1;

end


function [fullStep, c_sqr_sum, flag, c_hat, c_hat_l, c_hat_r, active, activity, n_active, min_ineq, m_x, m_x_try, m_rho, f_m, obj_val, lambda_t_c_hat] = line_search(fullStep, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, activity, n_active, c_hat, c_hat_l, c_hat_r, ineq_high, obj_val_old, param)
% adaptive Step size based on the status of convergence. 

m = 0;
no_activation = -1;
f_m = 0; f_l = 0; f_r = 0; c2s_l = 0; c2s_r = 0; bad_m = 0; bad_l = 0; bad_r = 0;
viol_l = 0; viol_r = 0; act_l = 0; act_m = 0; act_r = 0;
act_keep_l = no_activation;
act_keep_r = no_activation;
c_sqr_sum = 0;

[f_m, act_m, bad_m, c_sqr_sum, c_hat, min_ineq, m_x_try, obj_val] = aug_lagrangian(m, c_hat, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, param);
obj_diff = abs(obj_val_old - obj_val);
if obj_diff > 1e-2
    step_factor = 0.6;
elseif obj_diff > 1e-4
    step_factor = 0.5;
elseif obj_diff > 1e-5
    step_factor = 0.4;   
elseif obj_diff > 1e-6
    step_factor = 0.3;  
else
    step_factor = 0.2;  
end

stepSize = param.step_large;
for i = 1 : 100
    step = fullStep * stepSize;
    [f_l, act_l, bad_l, c2s_l, c_hat_l, min_ineq, m_x_try] = aug_lagrangian(m - step, c_hat_l, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, param);
    [f_r, act_r, bad_r, c2s_r, c_hat_r, min_ineq, m_x_try] = aug_lagrangian(m + step, c_hat_r, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, param);
    viol_l = act_l ~= no_activation || bad_l > param.constr_tol;
    viol_r = act_r ~= no_activation || bad_r > param.constr_tol;
    
    if (~viol_l && f_l < f_m && (viol_r || f_l < f_r))
        f_m = f_l; bad_m = bad_l; m = m - step; act_keep_r = no_activation;
    elseif (~viol_r && f_r < f_m && (viol_l || f_r < f_l))
        f_m = f_r; bad_m = bad_r; m = m + step; act_keep_l = no_activation;
    else
        act_keep_r = act_r; act_keep_l = act_l;
    end
    
    stepSize = stepSize * 0.5;
    stepSize = stepSize * step_factor;
    if stepSize <= param.line_tol
        break;
    end
end

[f_m, act_m, dum_bad, c_sqr_sum, c_hat, min_ineq, m_x_try, obj_val, lambda_t_c_hat] = aug_lagrangian(m, c_hat, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, param);
if f_m == 0
    flag_test = 1;
end

% m_x_try = m_x;
m_x = m_x_try;

if m~= 0.0
    m_rho = m_rho * (param.constr_tol * param.c1rho / (param.constr_tol * param.c0rho - bad_m) + param.c2rho) + param.rho_limit;
end

fprintf(" %6.0e %10.6f %8.1e %8.1e", bad_m, f_m, m_rho,  m);

[c_hat, active, activity, n_active, flag_l] = activate(act_keep_l, m_x_try, faces, c_hat, activity, active, n_active, "left");
[c_hat, active, activity, n_active, flag_r] = activate(act_keep_r, m_x_try, faces, c_hat, activity, active, n_active, "right");
% fprintf('\n%d %d %d %d\n', act_keep_l, act_keep_r, flag_l, flag_r);
if flag_l || flag_r
    flag = 1;
    return;
end

fullStep = fullStep * param.step_small;
if (abs(m) > fullStep)
    fullStep = abs(m);
end
flag = 0;

end


function [c_hat_try, constr_ineq, badness, min_ineq, m_x_try, act] = step_and_check(faces, step, m_x, m_dx, active, n_active, c_hat, ineq_high)

m_x_try = spher_step(step, m_x, m_dx);
[c_hat_try, constr_ineq] = constraints(m_x_try, faces);
c_hat_try = c_hat_try(1 : end - 1);
inequal_mat = reshape(constr_ineq', [], 1);

nface = size(faces, 1);
[c_hat_try, badness, min_ineq, act] = check_constraints(nface - 1, inequal_mat, c_hat_try, active, n_active, c_hat, ineq_high);

end


function [lagr, activate_, badness, c_sqr_sum, c_hat_try, min_ineq, m_x_try, obj_val, lambda_t_c_hat] = aug_lagrangian(step, c_hat, faces, neighbors, total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, param)
no_activation = -1;
nface = size(faces, 1);

[c_hat_try, ~, badness, min_ineq, m_x_try, activate_] = step_and_check(faces, step, m_x, m_dx, active, n_active, c_hat, ineq_high);
if(activate_ ~= no_activation)
    c_sqr_sum = 0;
    lagr = 0;
    lambda_t_c_hat = -m_lambda' * c_hat_try;
    obj_val = goal_func(m_x_try, neighbors, total_nb_num);
    return;
end

obj_val = goal_func(m_x_try, neighbors, total_nb_num);
lambda_t_c_hat = -m_lambda' * c_hat_try;
lagr = obj_val + lambda_t_c_hat;
c_hat_try_eq = c_hat_try(1 : nface - 1);
c_hat_try_ineq = c_hat_try(nface : end);
c_sqr_sum = sum(c_hat_try_eq .^ 2) + sum(c_hat_try_ineq(c_hat_try_ineq < 0) .^ 2);

lagr = lagr + m_rho * c_sqr_sum;
    
end


function [equal_mat, inequal_mat] = constraints(m_x, faces)

nface = size(faces, 1);
desired_area = 4 * pi / nface;

[equal_mat, inequal_mat] = spher_area4_function(m_x, faces);
equal_mat = equal_mat - desired_area;

end


function [inequal] = one_inequality(m_x, faces, act)

if nargin < 3
    act = 1 : size(faces, 1);
end

sine_id = rem(act + 4 - 1, 4) + 1;

face_id = ceil(act / 4);
[~, inequal_mat] = spher_area4_function(m_x, faces(face_id, :));
inequal = inequal_mat(sine_id);

end


function [c_hat_try, badness, min_ineq, act] = check_constraints(n_equal, inequal_mat, c_hat_try, active, n_active, c_hat, ineq_high)

min_ineq = 1;
min_pos = -1;

active = active(1 : end - 1);
c_hat_try(n_equal + 1 : n_equal + n_active) = inequal_mat(active);

inequal_mat(active) = min_ineq;
[min_val, min_ind] = min(inequal_mat);
if min_val < min_ineq
    min_ineq = min_val;
    min_pos = min_ind;
end

badness = 0;
worst = -1;

if min_ineq < -ineq_high
    act = min_pos;
    return;
end

increase_mat = c_hat_try .^ 2 - c_hat(1:n_equal + n_active) .^ 2;
[badness, worst] = max(increase_mat);
act = -1;

end


function [c_hat, active, activity, n_active] = activate_initial(c_hat, inequal_mat, ineq_high, no_activation)
% 03/12/2018 vectorize it 

sprintf('initially active: ');
n_inequal = numel(inequal_mat);
activity = zeros(n_inequal, 1);
inds = find(inequal_mat(:) < -ineq_high);
active = inds';
n_active = numel(active);
c_hat(end + 1 : end + n_active) = inequal_mat(inds);
activity(inds) = 6;

active(n_active + 1) = no_activation;
    
end 


function [m_x_1] = spher_step(step, m_x, m_dx)

m_x_1 = m_x + step .* m_dx;
m_x_1 = m_x_1 ./ sqrt(sum(m_x_1 .^ 2, 2));

end


function [areas, spat] = spher_area4_function(vertices, faces)

a = vertices(faces(:, 1), :);
b = vertices(faces(:, 2), :);
c = vertices(faces(:, 3), :);
d = vertices(faces(:, 4), :);

ab = dot(a, b, 2);
ac = dot(a, c, 2);
ad = dot(a, d, 2);
bc = dot(b, c, 2);
bd = dot(b, d, 2);
cd = dot(c, d, 2);

Ca = bd - ad .* ab;
Cb = ac - ab .* bc;
Cc = bd - bc .* cd;
Cd = ac - cd .* ad;

spat = zeros(size(faces, 1), 4);
spat(:, 1) = determinant_computing_3d(d, a, b);
spat(:, 2) = determinant_computing_3d(a, b, c);
spat(:, 3) = determinant_computing_3d(b, c, d);
spat(:, 4) = determinant_computing_3d(c, d, a);

areas = -atan2(Ca, spat(:, 1)) - atan2(Cb, spat(:, 2)) -atan2(Cc, spat(:, 3)) -atan2(Cd, spat(:, 4));
areas = mod(areas + 8.5 * pi, pi) - 0.5 * pi;

end


function [det_vec] = determinant_computing_3d(A, B, C)
% A, B, C are row vectors with dim n X 3

det_vec = A(:, 1) .* B(:, 2) .* C(:, 3) + A(:, 2) .* B(:, 3) .* C(:, 1) + A(:, 3) .* B(:, 1) .* C(:, 2) ...
          - A(:, 3) .* B(:, 2) .* C(:, 1) - A(:, 2) .* B(:, 1) .* C(:, 3) - A(:, 1) .* B(:, 3) .* C(:, 2);

end
