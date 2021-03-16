function [mix, options, errlog, mixlog] = ml_gmmem(mix, x, options, weights)
%GMMEM	EM algorithm for Gaussian mixture model.
%
%	Description
%	[MIX, OPTIONS, ERRLOG] = ML_GMMEM(MIX, X, OPTIONS) uses the Expectation
%	Maximization algorithm of Dempster et al. to estimate the parameters
%	of a Gaussian mixture model defined by a data structure MIX. The
%	matrix X represents the data whose expectation is maximized, with
%	each row corresponding to a vector.    The optional parameters have
%	the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%   OPTIONS(2) is for how to estimate the population with one component. If
%   it is 0, full covariance matrix will be used. If it is 1, the
%   covariance matrix will be the same as specified by MIX
%
%	OPTIONS(3) is a measure of the absolute precision required of the
%	error function at the solution. If the change in log likelihood
%	between two steps of the EM algorithm is less than this value, then
%	the function terminates.
%
%	OPTIONS(5) is set to 1 if a covariance matrix is reset to its
%	original value when any of its singular values are too small (less
%	than MIN_COVAR which has the value eps).   With the default value of
%	0 no action is taken.
%
%   OPTIONS(6) provides a minimum regularization value for any of the
%   diagonal values of the learned covariance matrix. Default value is
%   1e-5. If set to 0 no action is taken.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%
%	The optional return value OPTIONS contains the final error value
%	(i.e. data log likelihood) in OPTIONS(8).
%
%	See also
%	GMM, GMMINIT
%

%	Copyright (c) Ian T Nabney (1996-2001)
%   18-Jul-2006 Modified by T. Zhao

if ~exist('weights','var')
    weights = [];
end

if isempty(weights)
    weights = ones(size(x,1),1);
end

display = options(1);

%If there is only one component
if mix.ncentres==1
    errlog = 0;
    mixlog = {};
    mix.priors = 1;
    
    if options(2)==0
        [mix.centres,mix.covars] = ml_objgauss([x weights]);
        mix.covar_type = 'full';
    else
        [mix.centres,mix.covars] = ml_objgauss([x weights],mix.covar_type);
    end
    
    % Regularize covariance matrix
    mix = regularizeCovars(mix, options(6));
    
    return;
end

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);

if ~isempty(errstring)
    error(errstring);
end

[ndata, xdim] = size(x);
ndata = sum(weights);

% Sort out the options
if (options(14))
    niters = options(14);
else
    niters = 100;
end

store = 0;
if (nargout > 2)
    store = 1;	% Store the error values to return them
    errlog = zeros(1, niters);
end
if (nargout > 3)
    store = 2; % Store mixtures for each iteration
    mixlog = cell(1, niters);
end

test = 0;
if options(3) > 0.0
    test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
    if display >= 0
        disp('check_covars is on');
    end
    check_covars = 1;	% Ensure that covariances don't collapse
    MIN_COVAR = eps;	% Minimum singular value of covariance matrix
    init_covars = mix.covars;
end

intervaldata = 0;
if (options(15))
    intervaldata = 1;
end

scale = options(16);
if ~scale
    scale = 1;
end

es = [];
e = Inf;
% Main loop of algorithm
for n = 1:niters
    %added by Ivan E. Cao-Berg to handle non singular matrix error
    try
        % Calculate posteriors based on old parameters
        [post, act] = gmmpost(mix, x, scale);
    catch
        return
    end
    
    % Calculate error value if needed
    if (display | store | test)
        prob = act*(mix.priors)';
        % Error value is negative log likelihood of data
        e0 = e;
        e = - sum(log(prob).*weights);
        %         disp(['Step' int2str(n) ', Error value: ' num2str(e)])
        if abs((e0-e)/e) < 1e-6
            options(8) = -sum(log(gmmprob(mix, x)));
            return
        end
        %         es = [es e];
        %         if(length(es)>2)
        %             plot(es(3:end)-es(2:end-1));
        %             drawnow
        %         end
        if store
            errlog(n) = e;
        end
        if store > 1
            mixlog{n} = mix;
        end
        if display > 0
            fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
        end
        if test
            if (n > 1 & abs(e - eold) < options(3))
                options(8) = e;
                return;
            else
                eold = e;
            end
        end
    end
    
    % Adjust the new estimates for the parameters
    post = ml_multrow(post',weights')';
    new_pr = sum(post, 1);
    new_c = post' * x;
    
    % Now move new estimates to old parameter vectors
    mix.priors = new_pr ./ ndata;
    
    mix.centres = new_c ./ (new_pr' * ones(1, mix.nin));
    
    switch mix.covar_type
        case 'spherical'
            n2 = dist2(x, mix.centres);
            for j = 1:mix.ncentres
                v(j) = (post(:,j)'*n2(:,j));
            end
            clear n2
            mix.covars = ((v./new_pr))./mix.nin;
            % T. Peng +: Do not ignore digitalization
            if intervaldata
                mix.covars = mix.covars + 1/12*(1/scale)^2;
            end
            if check_covars
                % Ensure that no covariance is too small
                for j = 1:mix.ncentres
                    if mix.covars(j) < MIN_COVAR
                        mix.covars(j) = init_covars(j);
                    end
                end
            end
        case 'diag'
            for j = 1:mix.ncentres
                diffs = x - (ones(size(x,1), 1) * mix.centres(j,:));
                mix.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1, ...
                    mix.nin)), 1)./new_pr(j);
            end
            % T. Peng +: Do not ignore digitalization
            if intervaldata
                mix.covars = mix.covars + 1/12*(1/scale)^2;
            end
            if check_covars
                % Ensure that no covariance is too small
                for j = 1:mix.ncentres
                    if min(mix.covars(j,:)) < MIN_COVAR
                        mix.covars(j,:) = init_covars(j,:);
                    end
                end
            end
        case 'full'
            for j = 1:mix.ncentres
                diffs = x - (ones(size(x,1), 1) * mix.centres(j,:));
                diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
                mix.covars(:,:,j) = (diffs'*diffs)/new_pr(j);
            end
            % T. Peng +: Do not ignore digitalization
            if intervaldata
                mix.covars = mix.covars;% + 1/12*(1/scale)^2;
            end
            if check_covars
                % Ensure that no covariance is too small
                for j = 1:mix.ncentres
                    if min(svd(mix.covars(:,:,j))) < MIN_COVAR
                        mix.covars(:,:,j) = init_covars(:,:,j);
                    end
                end
            end
        case 'ppca'
            for j = 1:mix.ncentres
                diffs = x - (ones(size(x,1), 1) * mix.centres(j,:));
                diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
                [tempcovars, tempU, templambda] = ...
                    ppca((diffs'*diffs)/new_pr(j), mix.ppca_dim);
                if length(templambda) ~= mix.ppca_dim
                    error('Unable to extract enough components');
                else
                    mix.covars(j) = tempcovars;
                    mix.U(:, :, j) = tempU;
                    mix.lambda(j, :) = templambda;
                end
            end
            if check_covars
                if mix.covars(j) < MIN_COVAR
                    mix.covars(j) = init_covars(j);
                end
            end
        otherwise
            error(['Unknown covariance type ', mix.covar_type]);
    end
    
    clear post
    clear act
    
    mix = regularizeCovars(mix, options(6));
    
end

options(8) = -sum(log(gmmprob(mix, x)));
if (display >= 0)
    disp('Warning: Maximum number of iterations has been exceeded');
end

end

function mix = regularizeCovars(mix, regValue, debug)
% Regularizes the covariance matrices of mix to regValue

% Regularization value of 0 = don't do regularization
if (regValue == 0)
    return
end

% Don't do regularization for not full matrices
if (~strcmp(mix.covar_type, 'full'))
    return
end

if nargin < 3
    debug = false;
end

if ~checkForPositiveDefinite(mix)
    1;
end


d = size(mix.covars, 1);
for j = 1 : mix.ncentres
    
    if debug && (any(diag(mix.covars(:,:,j)) < regValue))
        disp(['Covariance #' int2str(j) 'below regularization value']);
    end
    
    % Regularize by setting any small eigenvalues to the regularization
    % value
    covar = mix.covars(:,:,j);
    [V, D] = eig(covar);
    D(logical(eye(d))) = max(regValue*ones(d,1), diag(D));
    mix.covars(:,:,j) = V*D/V;
end

if ~checkForPositiveDefinite(mix)
    error('Did not regularize adequately');
end

end



function allPosDef = checkForPositiveDefinite(mix, debug)
% Checks if covariances are all positive definite
% Inputs:
% mix: the learned mixture
% debug: (optional) whether to print if not all positive definite

allPosDef = true;

if nargin < 2
    debug = false;
end

% If square covariance matrices, check the full matrix
if size(mix.covars, 1) == size(mix.covars, 2)
    for j = 1:mix.ncentres
        [~,p] = chol(mix.covars(:,:,j));
        if p
            allPosDef = false;
        end
    end
% If linear covariance matrix, must be multicomponent spherical.
% Check that each variance is not 0.
elseif (size(mix.covars, 1) == 1)
    allPosDef = any(mix.covars < eps);
end


if debug
    if ~allPosDef
        disp('Not all learned covarinces are positive definite');
    end
end

end
