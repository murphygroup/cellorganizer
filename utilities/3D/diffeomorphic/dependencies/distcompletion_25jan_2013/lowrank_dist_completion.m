function [Y infos] = lowrank_dist_completion(algofun,Is,Js,knownDists,Y,params)
% [Y infos] = lowrank_dist_completion(algofun,Is,Js,knownDists,params)
%
% Perform low-rank distance completion 
% with automatic detection of the rank
%
% It implements the ideas of Journee et al. 2010 (SIMAX)
% applied to the problem of low-rank distance completion.
%
% See the paper Mishra et al. 2011 (CDC) for details
%
% Parameters
%
%   algofun      function handle on the optimization algorithm
%                (trust region or gradient descent)
%
%   Is           m-by-1 vector containing first indices of the points
%   Js           m-by-1 vector containing second indices of the points
%   knownDists   m-by-1 vector of known distances
%   Y            m-by-p0 initial condition of rank p0
%   params       structure array containing algorithm parameters
%                (see DefaultParams for details)
%
% Output
%   Y            n-by-r solution matrix of rank r
%   infos        structure array with additional information
%

% Authors:
% Bamdev Mishra and Gilles Meyer
% {b.mishra, g.meyer}@ulg.ac.be
% Last modified: 25 Jan, 2013.


    % Set default params
    params = DefaultParams(params);
    
    verb = params.verb;
    sig_A = params.sig_A;
    max_step = params.max_step;
    ls_maxiter = params.ls_maxiter;
    
    smin_tol = params.smin_tol;
    vp_tol = params.vp_tol;
    
    [n,p0] = size(Y);
    m = length(knownDists);
    
    if isfield(params,'pmax'),
        pmax = params.pmax;
    else
        pmax = n;
    end

    EIJ = speye(n);
    EIJ = EIJ(:,Is) - EIJ(:,Js);
    
    % Send this matrix to the inner algorithm 
    % so that it does not have to recompute it
    params.EIJ = EIJ;

    % Make eigs silent
    opts.disp = 0;    
    
    % Record cost function values and new rank locations
    infos.costs = [];
    infos.newRank = [];
    
    % print calling algorithm
    name_function = func2str(algofun);

    fprintf('>> Inner algorithm is %s <<\n',name_function);

    p = p0;
    while (p <= pmax), % When p=n a global min is attained for sure
        
        
            fprintf('>> Rank %d <<\n',p);
        
        
        if (p > p0),
            
            if isempty(restartDir), % If no restart dir avail. do random restart
                
                disp('No restart dir available, random restart is performed');
                Y = randn(n,p);
                
            else % Perform line-search based on the restart direction 
                
                disp('>> Line-search with restart direction');
                Y(:,p) = 0; % Append a column of zeroes
                
                Z = Y(Is,:) - Y(Js,:);
                estimDists = sum(Z.^2,2);
                errors = (estimDists - knownDists);
                                
%                 grad_Y = EIJ * sparse(1:m,1:m,2 * errors / m,m,m) * Z;
                
                costBefore = mean(errors.^2);
                fprintf('>> Cost before = %f\n',costBefore);

%                 step = max_step;
                step = abs(s_min);
                for i = 1:ls_maxiter,
                    
                    % Update
                    Y(:,p) = step*restartDir;
                    
                    % Compute cost
                    Z = Y(Is,:) - Y(Js,:);
                    estimDists = sum(Z.^2,2);
                    errors = (estimDists - knownDists);
                
                    costAfter = mean(errors.^2);
                    fprintf('>> Cost after = %f\n',costAfter);
                    
                    % Armijo condition
%                     armijo = (costAfter - costBefore) <= sig_A * step * (restartDir'*grad_Y(:,p));
                    armijo = (costAfter - costBefore) <= sig_A * step * abs(s_min);
                    if armijo,
                        break;
                    else
                        step = step/2;
                    end
                    
                end
                
                % Check for sufficient decrease
                if (costAfter >= costBefore) || abs(costAfter - costBefore) < 1e-8,
                    disp('Decrease is not sufficient, random restart');
                    Y = randn(n,p);
                end
                
            end    

        end
        
        % Run algorithm
        [Y infos_algo] = feval(algofun,Is,Js,knownDists,Y,params);

        infos.costs = [infos.costs;infos_algo.costs];
        infos.newRank = [infos.newRank;infos_algo.iter];
        
        % Evaluate gradient of the convex cost function (i.e. wrt X)
        Z = Y(Is,:) - Y(Js,:);
        estimDists = sum(Z.^2,2);
        errors = (estimDists - knownDists);

        % Dual variable
        Sy = EIJ * sparse(1:m,1:m,2 * errors / m,m,m) * EIJ';
        
        % Compute smallest algebraic eigenvalue of Sy,
        % this gives us a descent direction for the next rank (v)
        % as well as a way to control progress toward the global
        % optimum (s_min)
        [v, s_min] = eigs(Sy, 1, 'SA', opts);
        
        % To check whether Y is rank deficient
        vp = svd(Y);
        
        % Stopping criterion
        fprintf('>> smin = %f, and min(vp) = %f\n',s_min,min(vp));
        if (s_min  > -smin_tol) || (min(vp) < vp_tol),
            break;
        end
        
        p = p + 1;

        if (s_min < -1e-10),
            restartDir = v;
        else
            restartDir = [];
        end
        clear Sy v;

    end
    
    infos.newRank = cumsum(infos.newRank);

end
