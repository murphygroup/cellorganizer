function params = DefaultParams(params)
% Set default parameters for gradient descent and trust-region algorithms
  
% Authors:
% Bamdev Mishra and Gilles Meyer
% {b.mishra,g.meyer}@ulg.ac.be

  % Tolerance on absolute value
  if ~isfield(params,'tol'),
    params.tol = 1e-5;
  end
  
  % Tolerance on relative value
  if ~isfield(params,'vtol'),
    params.vtol = 1e-5;
  end
  
  % Tolerance on smallest eigenvalue of Sy, the dual variable
  if ~isfield(params,'smin_tol'),
    params.smin_tol = 1e-3;
  end
  
  % Tolerance on smallest singular value of Y
  if ~isfield(params,'vp_tol'),
    params.vp_tol = 1e-3;
  end  
  
  % Verbosity
  if ~isfield(params,'verb'),
    params.verb = true;
  end
  
  % Maximum step size during line-search
  if ~isfield(params,'max_step')
    params.max_step = 10;
  end

  % Parameter of the Armijo step
  if ~isfield(params,'sig_A'),
    params.sig_A = 0.5;
  end
  
  % Maximum number of line-search steps
  if ~isfield(params,'ls_maxiter'),
    params.ls_maxiter = 100;
  end
  
  % Maximum number of iterations
  if ~isfield(params,'maxiter'),
    params.maxiter = 10000;
  end
  
  % Maximum number of iterations for trust-region
  if ~isfield(params,'maxiter_tr'),
    params.maxiter_tr = 1000;
  end

end