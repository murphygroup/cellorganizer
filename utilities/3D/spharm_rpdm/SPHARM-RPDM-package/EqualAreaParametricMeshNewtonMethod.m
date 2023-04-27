function [m_x_out, cost_mat, is_success] = EqualAreaParametricMeshNewtonMethod(...
          vertices, faces, m_x, par)
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
% (Styner et al. 2006) which is based on Brechbuhler & Gerig et al. 1995.
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
% 05/03/2018 use stricter criterion for bast_largr_diff
% 09/21/2018 check matlab version to decide whether to use lsminnorm or not
%% these 4 changes got reverted by Khaled's change below
% 10/16/2020 R.F.Murphy disable MaxIter increase and add to verbose output
% 10/22/2020 R.F.Murphy initialize m_x_best_bu
% 12/4/2020 R.F. Murphy disable singular matrix warnings
% 2/1/2021 R.F. Murphy fix error when number of iter is less than 21
% 4/24/2023 R.F. Murphy fix missing <CR> on termination message
%%
% ------------------------------------------------------------------------
% Modifications June 2021 by Khaled Khairy: khaled.khairy@stjude.org
% St. Jude Children's Research Hospital
% Comments in code: Search for "KK"
% Specific modifications:
%
% - Speed up of Jacobian calculation ~ 5x
%       * Jacobian pattern is precalculated in a first step and then used
%         afterwards
%       * Call to calc_area_jacobian_sphere has been fully vectorized
%         Calls modified function "calc_area_jacobian_sphere_KK.m"
%
% - Speed up of linesearch calculation ~ 5x
%       * Precalculating quantities for neighbor counts 3-12
%       * Vectorization of for loop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 8/10/2022 R.F.Murphy reinstate changes from 10/16/2020 through 2/1/2021
%                   also set debug and verbose flags from par

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')

default_par.debug = false;
default_par.verbose = false;
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
default_par.lsq_bound = 1000; % bound of solutions of least squares in the update.

if nargin < 4
    par = struct();
end
par = process_options_structure(default_par, par);

verbose = par.verbose;

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

%% KK: use neighbors to pre-calculate quantities for linesearch

% % KK: generate set of indices into m_x shaped to vectorize sum
% % this is possible because the mesh topology is constant
% % Since we have either three four five or six vertices we will have to do this
% % separately
%
% % determine neighbors{ix} number of vertices
nvert = size(neighbors, 1);
indx_rec.indx12 = zeros(nvert,1, 'logical');
indx_rec.indx11 = zeros(nvert,1, 'logical');
indx_rec.indx10 = zeros(nvert,1, 'logical');
indx_rec.indx9 = zeros(nvert,1, 'logical');
indx_rec.indx8 = zeros(nvert,1, 'logical');
indx_rec.indx7 = zeros(nvert,1, 'logical');
indx_rec.indx6 = zeros(nvert,1, 'logical');
indx_rec.indx5 = zeros(nvert,1, 'logical');
indx_rec.indx4 = zeros(nvert,1, 'logical');
indx_rec.indx3 = zeros(nvert,1, 'logical');
%numindx = zeros(nvert,1);
for ix = 1:nvert
    numindx(ix) = numel(neighbors{ix});    %disp(numindx(ix);
    indx_rec.indx12(ix) = numel(neighbors{ix})==12;
    indx_rec.indx11(ix) = numel(neighbors{ix})==11;
    indx_rec.indx10(ix) = numel(neighbors{ix})==10;
    indx_rec.indx9(ix) = numel(neighbors{ix})==9;
    indx_rec.indx8(ix) = numel(neighbors{ix})==8;
    indx_rec.indx7(ix) = numel(neighbors{ix})==7;
    indx_rec.indx6(ix) = numel(neighbors{ix})==6;
    indx_rec.indx5(ix) = numel(neighbors{ix})==5;
    indx_rec.indx4(ix) = numel(neighbors{ix})==4;
    indx_rec.indx3(ix) = numel(neighbors{ix})==3;
end
if max(numindx)>12, error('Maximum number of edges per vertex exceeded');end
indx_rec.mxind12 = zeros(12,3, sum(indx_rec.indx12)); % store indices for "six" vertices
indx_rec.mxind11 = zeros(11,3, sum(indx_rec.indx11)); % store indices for "five" vertices
indx_rec.mxind10 = zeros(10,3, sum(indx_rec.indx10)); % store indices for "six" vertices
indx_rec.mxind9 = zeros(9,3, sum(indx_rec.indx9)); % store indices for "five" vertices
indx_rec.mxind8 = zeros(8,3, sum(indx_rec.indx8)); % store indices for "four" vertices
indx_rec.mxind7 = zeros(7,3, sum(indx_rec.indx7)); % store indices for "three" vertices
indx_rec.mxind6 = zeros(6,3, sum(indx_rec.indx6)); % store indices for "six" vertices
indx_rec.mxind5 = zeros(5,3, sum(indx_rec.indx5)); % store indices for "five" vertices
indx_rec.mxind4 = zeros(4,3, sum(indx_rec.indx4)); % store indices for "four" vertices
indx_rec.mxind3 = zeros(3,3, sum(indx_rec.indx3)); % store indices for "three" vertices

count12 = 1;
count11 = 1;
count10 = 1;
count9 = 1;
count8 = 1;
count7 = 1;
count6 = 1;
count5 = 1;
count4 = 1;
count3 = 1;
for ix = 1:nvert
    if indx_rec.indx3(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [3 3]);
        indx_rec.mxind3(:,:,count3) = n_indmx;
        count3 = count3 + 1;
    elseif indx_rec.indx4(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [4 3]);
        indx_rec.mxind4(:,:,count4) = n_indmx;
        count4 = count4 + 1;
    elseif indx_rec.indx5(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [5 3]);
        indx_rec.mxind5(:,:,count5) = n_indmx;
        count5 = count5 + 1;
    elseif indx_rec.indx6(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [6 3]);
        indx_rec.mxind6(:,:,count6) = n_indmx;
        count6 = count6 + 1;
    elseif indx_rec.indx7(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [7 3]);
        indx_rec.mxind7(:,:,count7) = n_indmx;
        count7 = count7 + 1;
    elseif indx_rec.indx8(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [8 3]);
        indx_rec.mxind8(:,:,count8) = n_indmx;
        count8 = count8 + 1;
    elseif indx_rec.indx9(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [9 3]);
        indx_rec.mxind9(:,:,count9) = n_indmx;
        count9 = count9 + 1;
    elseif indx_rec.indx10(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [10 3]);
        indx_rec.mxind10(:,:,count10) = n_indmx;
        count10 = count10 + 1;
    elseif indx_rec.indx11(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [11 3]);
        indx_rec.mxind11(:,:,count11) = n_indmx;
        count11 = count11 + 1;
    elseif indx_rec.indx12(ix)
        n_indmx = reshape([neighbors{ix}'; neighbors{ix}'+nvert; neighbors{ix}' + 2*nvert], [12 3]);
        indx_rec.mxind12(:,:,count12) = n_indmx;
        count12 = count12 + 1;
    end
end
indx_rec.nb_order = [find(indx_rec.indx3); find(indx_rec.indx4);...
                     find(indx_rec.indx5); find(indx_rec.indx6);...
                     find(indx_rec.indx7); find(indx_rec.indx8);...
                     find(indx_rec.indx9); find(indx_rec.indx10);...
                     find(indx_rec.indx11); find(indx_rec.indx12)];
% % mxind3, mxind4, mxind5 and mxind6 (and so on) store indices into m_x
% % to get the sum we can use for example: xx = [squeeze([sum(m_x(mxind3),1)])]';
%%

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
% 8/10/2022 reinstate commenting out of this increase
% may change MaxIter if nvert is very large
%MaxIter = max(MaxIter, ceil(nvert / 50));
best_largr = 100000;
largr_last = best_largr;
best_largr_bu = best_largr;
best_largr_diff = 10;
m_x_best = m_x;
m_x_best_bu = m_x;
cost_mat = zeros(MaxIter, 3);
obj_val = 1000000;
lambda_t_c_hat = 1e5;
% check whether the jacobi matrix ill
is_ill_ever = true;
last_ill_iter = 1;
ill_dist = 5;

% set up optimization options for lsqlin
% optm_options = optimoptions('lsqlin','Algorithm', 'trust-region-reflective', 'Display','none', 'FunctionTolerance', 1e-8);

% KK: Call once to generate index matrix cmx_ind
[jacobi_mat, is_jacobi_ill, cmx_ind] = calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, c_hat, J_sines_func, par);


% iterate to update
large = 900000000;

%check plot flag for KK
if ~isfield(par, 'plot_flag_nm')
    par.plot_flag_nm = 0;
end
if verbose
    fprintf("In EqualAreaParametricMeshNewtonMethod, MaxIter=%f\n",MaxIter);
    fprintf("Iter,nwtnstep,actkeep,gamma,ngradZ,csqrsum,cost,ineghigh,largrdiff\n");
end
for iter = 1 : MaxIter
    if par.plot_flag_nm == 1
        %% KK: visualize parameterization progress
        clf;set(gcf,'renderer','opengl','color','k', 'InvertHardCopy', 'off');
        patch('Vertices', m_x, 'Faces', faces,'FaceColor', 'b', 'EdgeColor','k');
        %mv = surface_mesh(m_x, faces);plot(mv);
        axis equal;axis off;
        view(60, 38);axis tight;
        title([num2str(iter) ': Parameterization optimization run'], 'Color', 'w');
        drawnow;
    end
    m_x_out =  m_x;
    if verbose, fprintf(" %3d", iter);   end
    
    % jacobi_mat_1 = calc_jacob_matrix_func_1(vertices, faces, neighbors, face_ids, m_x, c_hat, active, par);
    [jacobi_mat, is_jacobi_ill] = calc_jacob_matrix_func_2(vertices, faces,...
        m_x, active, n_active, c_hat, J_sines_func, par, cmx_ind);
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
        if verLessThan('matlab', '9.4')
            m_newton_dir = box_constrained_lsq_admm(jacobi_mat, c_hat, lsq_bound, -lsq_bound);
        else
            m_newton_dir = lsqminnorm(jacobi_mat, c_hat);
            if any(abs(m_newton_dir) > lsq_bound)
                m_newton_dir = box_constrained_lsq_admm(jacobi_mat, c_hat, lsq_bound, -lsq_bound);
            end
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
    
    if verbose, fprintf(" %6.0e %6d", newtonStep, act_keep);end
    m_x = spher_step(newtonStep, m_x, m_dx);
    c_hat_temp = c_hat;
    c_hat = c_hat_l;
    c_hat_l = c_hat_temp;
    c_hat_r = c_hat;
    
    grad = calc_gradient_func(m_x, neighbors);
    grad = reshape(grad', [], 1);
    
    [jacobi_mat, is_jacobi_ill] = calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, c_hat, J_sines_func, par, cmx_ind);
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
        if verLessThan('matlab', '9.4')
            m_lambda = box_constrained_lsq_admm(jacobi_mat', grad, lsq_bound, -lsq_bound);
        else
            m_lambda =  lsqminnorm(jacobi_mat', grad);
            if any(abs(m_lambda) > lsq_bound)
                m_lambda = box_constrained_lsq_admm(jacobi_mat', grad, lsq_bound, -lsq_bound);
            end
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
            [c_hat, flag, activity, active, n_active] = inactivate(nface, c_hat, activity, active, n_active, i, verbose);
            inactivated = inactivated + flag;
        end
    end
    
    if inactivated > 0
        [c_hat, active, activity, n_active, flag] = activate(act_keep, m_x_try, faces, c_hat, activity, active, n_active, "in_out", verbose);
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
    if verbose, fprintf(" %6.3f %8.1e", gamma, sqrt(norm2gradZ));end
    
    [m_alpha_step, c_sqr_sum, activated, c_hat, c_hat_l, c_hat_r, active,...
        activity, n_active, min_ineq, m_x, m_x_try, m_rho, largr, obj_val,...
        lambda_t_c_hat] = ...
        line_search(m_alpha_step, faces, neighbors, total_nb_num, m_x, ...
        m_dx, m_lambda, m_rho, active, activity, n_active, c_hat, c_hat_l,...
        c_hat_r, ineq_high, obj_val, par, indx_rec);
    
    [c_hat, active, activity, n_active, flag] = ...
        activate(act_keep, m_x_try, faces, c_hat, activity, active,...
        n_active, "newton", verbose);
    
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
    if verbose, fprintf(" %10.3e %10.3e %10.3e %10.3e\n", sqrt(c_sqr_sum), cost, ineq_high, largr_diff);end
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
    % 06/25/2018 add check of previous iteration for cost
    if iter > 2
        if (...
        (cost_mat(iter -1, 3) < (10 * par.cost_tol) && cost < par.cost_tol)...
        ||...
        (largr > 0 && ...
        abs(cost_mat(iter-1, 1) - cost_mat(iter-2, 1)) < (10 * par.largr_tol) &&...
        abs(largr_diff) < par.largr_tol)...
        )
            cost_mat(iter + 1 : end, :) = [];
            m_x_best = m_x;
            if verbose fprintf("\n"); end
            fprintf("Exiting with acceptable cost after %d iterations\n",iter);
            break;
        end
    end
    %     toc;
    
    if par.debug && rem(iter, 1000) == 0
        save('workspace.mat')
    end
end

if cost > 1e-2 || abs(lambda_t_c_hat) >  1e-2
    is_success = false;
else
    is_success = true;
end
if abs(best_largr_diff) > 1e-3
    is_success = false;
    m_x = m_x_best_bu;
else
    m_x = m_x_best;
end
% 05/05/2018 use largargian seq as a criterion
%if sum(abs(cost_mat(end-20:end, 1) - cost_mat(end, 1)) < 1e-4) < 14
% check if enough iterations gave small enough change in cost_mat
if size(cost_mat,1) > 20
    last20 = abs(cost_mat(end-20:end,1) - cost_mat(end,1));
    last20notnan = find(~isnan(last20));
    if sum(last20(last20notnan)<1e-4) < 14
        is_success = false;
        m_x = m_x_best_bu;
    else
        m_x = m_x_best;
    end
    m_x = m_x_best;
end


if verbose, disp('Optimization complete!'); end

if par.debug
    save('workspace.mat')
end

warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')

end


function [obj] = goal_func(m_x, neighbors, total_nb_num, ...
    indx_rec)

% KK: This is a slow loop. Now uses precalculated quantities
% nvert = size(m_x, 1);
% nbsum_vec = zeros(nvert, 3);
% for i = 1 : nvert
%     nbsum_vec(i, :) = sum(m_x(neighbors{i}, :));
% end

%%% KK: generate vectorized nbsum_vec and restore original row order of nbsum_vec

v3 = squeeze(sum(m_x(indx_rec.mxind3),1))';if numel(v3)==3, v3=v3';end
v4 = squeeze(sum(m_x(indx_rec.mxind4),1))';if numel(v4)==3, v4=v4';end
v5 = squeeze(sum(m_x(indx_rec.mxind5),1))';if numel(v5)==3, v5=v5';end
v6 = squeeze(sum(m_x(indx_rec.mxind6),1))';if numel(v6)==3, v6=v6';end
v7 = squeeze(sum(m_x(indx_rec.mxind7),1))';if numel(v7)==3, v7=v7';end
v8 = squeeze(sum(m_x(indx_rec.mxind8),1))';if numel(v8)==3, v8=v8';end
v9 = squeeze(sum(m_x(indx_rec.mxind9),1))';if numel(v9)==3, v9=v9';end
v10 = squeeze(sum(m_x(indx_rec.mxind10),1))';if numel(v10)==3, v10=v10';end
v11 = squeeze(sum(m_x(indx_rec.mxind11),1))';if numel(v11)==3, v11=v11';end
v12 = squeeze(sum(m_x(indx_rec.mxind12),1))';if numel(v12)==3, v12=v12';end

V = [v3; v4;v5;v6;v7; v8; v9;v10; v11; v12];

nbsum_vec = [indx_rec.nb_order V];
                   
nbsum_vec = sortrows(nbsum_vec);
nbsum_vec = nbsum_vec(:, 2:4);


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


function [A, is_A_ill, cmx_indx] = ...
    calc_jacob_matrix_func_2(vertices, faces, m_x, active, n_active, ...
    c_hat, J_sines_func, param, cmx_indx)
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

%%% KK: inspect
% m = surface_mesh(m_x, faces);plot(m);
%% KK: Generate a jacobian pattern and parallelize or vectorize the Jacobian calculation in the loop

% KK: get indices into m_x by preparing cmx_indx if not provided
if nargin<9 || isempty(cmx_indx)
    cmx_indx = zeros(3,4,nface);  % store indices into m_x
    for ix = 1 : nface
        cur_face = faces(ix, :);
        cmx_indx(:,:,ix) = (cur_face - 1) * 3 + (1:3)';
    end
end

% KK: generate an updated cmx for current configuration
cmx = m_x(cmx_indx);
% KK: calculate J_area vectorized
J_area_vec = calc_area_jacobian_sphere_KK(cmx);

% KK: Loop over faces as before, but without area jacobian calculation each
% time

for i = 1 : nface
    cur_face = faces(i, :);
    col_elm_ind = (cur_face - 1) * 3 + (1:3)';
    col_elm_ind = col_elm_ind(:);
    cmx = m_x(:, cur_face);
    
    % KK: Prepared J_area_vec by vectrization of calc_area_jacobian_sphere
    J_area = J_area_vec(:,:,i);
    
    %%%% KK: Very slow to call once for every face for every iteration ---
    %%%% use vectorized version instead of below
    % J_area = calc_area_jacobian_sphere(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4));
    
    
    
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
%%
% put pairs values to the matrix
index_val_mat = index_val_mat(1 : counter-1, :);
A = sparse(index_val_mat(:, 1), index_val_mat(:, 2), index_val_mat(:, 3), nface + n_active,  nvert * 3);
A(nface, :) = [];

end


function [c_hat, active, activity, n_active, flag] = activate(act, m_x_try, faces, c_hat, activity, active, n_active, info, verbose)
no_activation = -1;
flag = 0;
if act == no_activation
    return;
end
if verbose
    fprintf("\n]]] %s : constraint %d -> level %d ", info, act, activity(act) + 1);
end
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

if verbose
    fprintf(">>> %s activates constraint %d (pos %d); now active: %d\n                                                 ",  ...
    info, act, i, n_active);
end

flag = 1;

end


function [c_hat, flag, activity, active, n_active] = inactivate(nface, c_hat, activity, active, n_active, pos, verbose)

flag = 0;
if verbose
    fprintf('\n[[[ constraint %d -> level %d ;          ', active(pos), activity(active(pos)) - 1);
end

activity(active(pos)) = activity(active(pos)) - 1;
if activity(active(pos)) ~= 3
    return;
end

activity(active(pos)) = 2;
n_active = n_active - 1;
if verbose
    fprintf('<<< inactivate constraint %d at position %d; now active: %d\n', active(pos), pos, n_active);
end

c_hat(nface - 1 + pos) = [];
active(pos) = [];
flag = 1;

end


function [fullStep, c_sqr_sum, flag, c_hat, c_hat_l, c_hat_r, active, ...
    activity, n_active, min_ineq, m_x, m_x_try, m_rho, f_m, obj_val, ...
    lambda_t_c_hat] = ...
    line_search(fullStep, faces, neighbors, total_nb_num, m_x, m_dx, ...
    m_lambda, m_rho, active, activity, n_active, c_hat, c_hat_l, c_hat_r,...
    ineq_high, obj_val_old, param,indx_rec)
% adaptive Step size based on the status of convergence.
verbose = false;
m = 0;
no_activation = -1;
f_m = 0; f_l = 0; f_r = 0; c2s_l = 0; c2s_r = 0; bad_m = 0; bad_l = 0; bad_r = 0;
viol_l = 0; viol_r = 0; act_l = 0; act_m = 0; act_r = 0;
act_keep_l = no_activation;
act_keep_r = no_activation;
c_sqr_sum = 0;


% KK: aug_lagrangian is a limiting step
[f_m, act_m, bad_m, c_sqr_sum, c_hat, min_ineq, m_x_try, obj_val] = ...
    aug_lagrangian(m, c_hat, faces, neighbors, total_nb_num, m_x, m_dx, ...
    m_lambda, m_rho, active, n_active, ineq_high, param, indx_rec);


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
    [f_l, act_l, bad_l, c2s_l, c_hat_l, min_ineq, m_x_try] = ...
        aug_lagrangian(m - step, c_hat_l, faces, neighbors, total_nb_num,...
        m_x, m_dx, m_lambda, m_rho, active, n_active, ineq_high, ...
        param, indx_rec);
    
    [f_r, act_r, bad_r, c2s_r, c_hat_r, min_ineq, m_x_try] = ...
        aug_lagrangian(m + step, c_hat_r, faces, neighbors, ...
        total_nb_num, m_x, m_dx, m_lambda, m_rho, active, n_active,...
        ineq_high, param, indx_rec);
    viol_l = act_l ~= no_activation || bad_l > param.constr_tol;
    viol_r = act_r ~= no_activation || bad_r > param.constr_tol;
    
    if (~viol_l && f_l < f_m && (viol_r || f_l < f_r))
        f_m = f_l; bad_m = bad_l; m = m - step; act_keep_r = no_activation;
    elseif (~viol_r && f_r < f_m && (viol_l || f_r < f_l))
        f_m = f_r; bad_m = bad_r; m = m + step; act_keep_l = no_activation;
    else
        act_keep_r = act_r; act_keep_l = act_l;
    end
    
    % stepSize = stepSize * 0.5;
    stepSize = stepSize * step_factor;
    if stepSize <= param.line_tol
        break;
    end
end

[f_m, act_m, dum_bad, c_sqr_sum, c_hat, min_ineq, m_x_try, obj_val,...
    lambda_t_c_hat] = ...
    aug_lagrangian(m, c_hat, faces, neighbors, ...
    total_nb_num, m_x, m_dx, m_lambda, m_rho, active, ...
    n_active, ineq_high, param, indx_rec);
if f_m == 0
    flag_test = 1;
end

% m_x_try = m_x;
m_x = m_x_try;

if m~= 0.0
    m_rho = m_rho * (param.constr_tol * param.c1rho / (param.constr_tol * param.c0rho - bad_m) + param.c2rho) + param.rho_limit;
end

if verbose, fprintf(" %6.0e %10.6f %8.1e %8.1e", bad_m, f_m, m_rho,  m);end

[c_hat, active, activity, n_active, flag_l] = activate(act_keep_l, m_x_try, faces, c_hat, activity, active, n_active, "left", verbose);
[c_hat, active, activity, n_active, flag_r] = activate(act_keep_r, m_x_try, faces, c_hat, activity, active, n_active, "right", verbose);
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


function [c_hat_try, constr_ineq, badness, min_ineq, m_x_try, act] = ...
    step_and_check(faces, step, m_x, m_dx, active, n_active, c_hat, ineq_high)

m_x_try = spher_step(step, m_x, m_dx);
[c_hat_try, constr_ineq] = constraints(m_x_try, faces);
c_hat_try = c_hat_try(1 : end - 1);
inequal_mat = reshape(constr_ineq', [], 1);

nface = size(faces, 1);
[c_hat_try, badness, min_ineq, act] = check_constraints(nface - 1, inequal_mat, c_hat_try, active, n_active, c_hat, ineq_high);

end


function [lagr, activate_, badness, c_sqr_sum, c_hat_try, min_ineq,...
    m_x_try, obj_val, lambda_t_c_hat] = ...
    aug_lagrangian(step, c_hat, faces, neighbors, total_nb_num, m_x, m_dx,...
    m_lambda, m_rho, active, n_active, ineq_high, ...
    param, indx_rec)
no_activation = -1;
nface = size(faces, 1);

[c_hat_try, ~, badness, min_ineq, m_x_try, activate_] = ...
    step_and_check(faces, step, m_x, m_dx, active, n_active, ...
    c_hat, ineq_high);
if(activate_ ~= no_activation)
    c_sqr_sum = 0;
    lagr = 0;
    lambda_t_c_hat = -m_lambda' * c_hat_try;
    obj_val = goal_func(m_x_try, neighbors, total_nb_num, ...
        indx_rec);
    return;
end

obj_val = goal_func(m_x_try, neighbors, total_nb_num, indx_rec);

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
