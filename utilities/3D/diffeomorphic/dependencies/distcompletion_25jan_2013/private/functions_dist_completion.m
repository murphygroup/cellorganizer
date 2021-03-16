function varargout=functions_dist_completion(varargin)
% Set of functions that are specific to distance completion problems
%
% Usages:
%
% Cost function:        f = functions_dist_completion('f',x,param);
% Euclidean Gradient:   grad = functions_dist_completion('grad_f',x,param);
% Euclidean Hessian:    hess = functions_dist_completion('hessian',x,eta,param);
%
% x is the current iterate
% eta is a search direction
% param is a cell array that contains {I,J,trueDists,EIJ} (see lowrank_dist_completion.m)

% Authors:
% Bamdev Mishra and Gilles Meyer
% {b.mishra,g.meyer}@ulg.ac.be

  fun_name = varargin{1};

  switch (fun_name),   

    % Objective
    case 'f',                   

      x=varargin{2};
      param=varargin{3};

      trueDists = param{3};
      EIJ = param{4};

      xij = EIJ'*x;

      estimDists = sum(xij.^2,2);

      varargout{1} = mean((estimDists - trueDists).^2);      

    % gradient on the Euclidean space R^{n \times p}
    case 'grad_f',                 

      x=varargin{2};
      param=varargin{3};
      
      trueDists = param{3};
      EIJ = param{4};

      m = length(trueDists);

      xij = EIJ'*x;
      
      estimDists = sum(xij.^2,2);

      varargout{1} = EIJ * sparse(1:m,1:m,2 * (estimDists - trueDists) / m, m, m) * xij;

    % Hessian in the direction eta on the Euclidean space R^{n \times p}
    case 'hessian',             

      x=varargin{2};
      eta=varargin{3};
      param=varargin{4};

      trueDists = param{3};
      EIJ = param{4};

      m = length(trueDists);

      xij = EIJ'*x;
      zij = EIJ'*eta;

      estimDists = sum(xij.^2,2);
      crossYZ = 2*sum(xij .* zij,2);

      varargout{1} = (EIJ*sparse(1:m,1:m,2 * (estimDists-trueDists) / m,m,m))*zij + (EIJ*sparse(1:m,1:m,2 * crossYZ / m,m,m))*xij;

  end

end