function [Y infos] = tr_dist_completion(Is,Js,knownDists,Y,params)
% [Y infos] = tr_dist_completion(Is,Js,knownDists,Y,params)
%
% Given (i,j,d)'s contained in [Is,Js,knownDists], the function compute a
% scalar product matrix X = Y*Y' that minimizes the mean quadratic error on
% observed distances. i and j's value range from 1 to n, the total number of points. 
%
% The optimization algorithm is trust-region.
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

  params = DefaultParams(params); % collect the parameters

  fun_obj=@functions_dist_completion;
  fun_set=@functions_manifold;

  param_tr{1} = Is;
  param_tr{2} = Js;
  param_tr{3} = knownDists;

  EIJ = speye(size(Y,1));
  param_tr{4} = EIJ(:,Is) - EIJ(:,Js);

  [Y, f, infos.costs] = trust_region(fun_set,fun_obj,Y,param_tr,params.tol,params.vtol,params.maxiter_tr,params.verb);
  infos.iter = length(infos.costs);

end
