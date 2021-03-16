function varargout=functions_manifold(varargin)
% Set of functions that defines geometrical objects on the quotient manifold 
% R^{n \times p}/O(p) where O_p is the orthogonal group
%
% Usages:
%
% Metric:               g = functions_manifold('metric',eta,zeta);
% Projection on H:      eta_h = functions_manifold('project',x,eta);
% Retraction:           x+ = functions_manifold('retraction',x,eta);
% Cost function:        f = functions_manifold('f',x,@fun_problem,param);
% Riemannian Gradient:  grad = functions_manifold('grad_f',x,@fun_problem,param);
% Riemannian Hessian:   hess = functions_manifold('hessian',x,@fun_problem,param);
%
% x is the current iterate
% eta is a search direction
% @fun_problem is a function handler on problem specific functions (see functions_dist_completion.m)
% param is a cell array that contains {I,J,trueDists,EIJ}

% Authors:
% Bamdev Mishra and Gilles Meyer
% {b.mishra,g.meyer}@ulg.ac.be


type=varargin{1};
fun_set=@functions_manifold;

switch type,
  
    case 'metric',              
        eta=varargin{3};
        zeta=varargin{4};
        varargout{1} = trace(eta'*zeta);

    case 'project',
      
        x=varargin{2};
        eta=varargin{3};
      
        SS=x'*x;           
        mat=x'*eta;
        AS=mat-mat';  
        O = sylvester2(SS,AS);

        varargout{1} = eta-x*O;
    
    case 'retraction',         
        
        x=varargin{2};
        eta=varargin{3};
        
        varargout{1} = x + eta;
        
    case 'f',                   

        x=varargin{2};
        fun_obj=varargin{3};
        param=varargin{4};        

        varargout{1} = feval(fun_obj,'f',x,param);              
        
    case 'grad_f',      

        x=varargin{2};
        fun_obj=varargin{3};
        param=varargin{4};

        varargout{1} = feval(fun_obj,'grad_f',x,param);

    case 'hessian',
        
        x=varargin{2};
        eta=varargin{3};
        fun_obj=varargin{4};
        param=varargin{5};   

        hess=feval(fun_obj,'hessian',x,eta,param);
        
        varargout{1}=feval(fun_set,'project',x,hess);
    
end

