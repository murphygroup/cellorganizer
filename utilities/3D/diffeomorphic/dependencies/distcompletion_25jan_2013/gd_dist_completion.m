function [Y infos] = gd_dist_completion(Is,Js,knownDists,Y,params)
% [Y infos] = gd_dist_completion(Is,Js,knownDists,Y,params)
%
% Given (i,j,d)'s contained in [Is,Js,knownDists], the function compute a
% scalar product matrix X = Y*Y' that minimizes the mean quadratic error on
% observed distances. i and j's value range from 1 to n, the total number
% of points.
%
% The optimization algorithm is gradient descent using the Armijo
% rule to select automatically the stepSize size.
%
%
% Parameters
%
%   Is           m-by-1 vector containing first indices of the points
%   Js           m-by-1 vector containing second indices of the points
%   knownDists   m-by-1 vector of known distances
%   Y            n-by-p initial condition matrix of rank p
%   params       structure array containing algorithm parameters
%                (see DefaultParams for details)
%
% Output
%   Y            n-by-p final matrix of rank p
%   infos        structure array with additional information
%

% Authors:
% Bamdev Mishra and Gilles Meyer
% {b.mishra,g.meyer}@ulg.ac.be
%
% Last modified: Jan 25, 2013


% Problem size
n = size(Y,1);

% Data size
m = length(knownDists);

% Collect parameters
params = DefaultParams(params);

% Shortcut parameter names
maxiter = params.maxiter;
sig_A = params.sig_A;
ls_maxiter = params.ls_maxiter;
max_step = params.max_step;
tol = params.tol;
vtol = params.vtol;
verb = params.verb;

% Fixed matrix that results from the gradient.
% For efficiency, it is computed once and for all at the beginning.
if isfield(params,'EIJ'),
    EIJ = params.EIJ;
else
    EIJ = speye(n);
    EIJ = EIJ(:,Is) - EIJ(:,Js);
end

% Compute initial cost
Z = EIJ'*Y;
estimDists = sum(Z.^2,2);
errors = (estimDists - knownDists);

cost = mean(errors.^2);

if verb,
    fprintf('[%0.4d] cost = %g\n',0,cost);
end

infos.costs = zeros(maxiter+1,1);
infos.costs(1) = cost;
infos.sqnormgrad = [];
infos.linesearch = [];



% Loop through all iterations
for iter = 1 : maxiter,
    
    S = sparse(1:m,1:m,2 * errors / m,m,m,m);
    
    % Compute gradient
    grad_Y = EIJ * S * Z;
    
    gradientNorm = norm(grad_Y,'fro')^2;
    infos.sqnormgrad = [infos.sqnormgrad ; gradientNorm];
    
    if iter > 1,
%         % Nocedal & Wright stepsize search
%         % A factor of 2 is used to avoid slow covergence due to underestimation of stepsize 
%         max_step = 2*stepSize * infos.sqnormgrad(end - 1) / gradientNorm;
        
        % Adaptive stepsize search as proposed by Mishra et al., arXiv:1209.0430
        % On an average we intend to perform only 2 linesearch per iteration.
        if s == 1,
            max_step = 2*max_step;
        elseif s >= 3,
            max_step = 2*stepSize;
        else % here s == 2,
            max_step = max_step; % Do nothing 
        end
        
    end
    
    % Perform line-search
    Yt = Y;
    stepSize = max_step;

    for s = 1 : ls_maxiter,
        
        % Update parameter
        Y = Yt - stepSize * grad_Y;
        
        % Evaluate new cost
        Z = EIJ'*Y;
        estimDists = sum(Z.^2,2);
        errors = (estimDists - knownDists);
        
        newcost = mean(errors.^2);
        
        % Check Armijo condition
        armijo = (cost - newcost) >= sig_A * stepSize * gradientNorm;
        if armijo,
            break;
        else
            stepSize = stepSize / 2;
        end
        
    end

    
    infos.costs(iter+1) = newcost;
    infos.linesearch= [infos.linesearch; s];
    if verb,
        fprintf('[%0.4d] cost = %g, grad norm = %g, #linesearch = %g\n',iter,newcost,gradientNorm, s);
    end
    
    % Stopping criterion
    if newcost < tol || (cost - newcost)/cost < vtol,
        break;
    end
    
    cost = newcost;
    
end


infos.costs = infos.costs(1:iter+1);
infos.iter = iter;

if iter >= maxiter,
    warning('MATLAB:MaxIterReached','Maximum number of iterations has been reached.');
end

end
