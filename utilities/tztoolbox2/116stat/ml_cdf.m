function y = ml_pdf(x,param)
%ML_PDF Computes a probablity density or mass function.
%   Y = ML_CDF(X,PARAM) returns the cumulative probablity or probablity 
%   density of the data X. It only supports univariate distributions.
%   PARAM is a structure ([pdf]) with the following fields: 
%       'name' - the name of the probablity density or mass function. It
%           can be:
%           'exp' : exponential function. 
%               function: f(x) = beta*exp(-beta*x), univariate
%               'beta' - location parameter.
%           'norm' : normal distribution.
%               x ~ N(mu,sigma), univariate
%               'mu' - mean
%               'sigma' - the standard deviation
%           'poiss' : poisson distribution.
%               function: f(x) = exp(-lamda)*lamda^x/x!, univariate, pmf
%               'lamda' - the shape parameter
%           'gamma' : gamma distribution
%               function: 
%                  f(x) = x^(alpha-1)*(beta^alpha*exp(-beta*x)/gamma(alpha)
%               'alpha' - shape parameter
%               'beta' - rate parameter
%           'mix' : Mixture distribution
%               f(x) = sum_i(p_i*g(x;param_i)), univaraite
%               'ps' - prior probabilities, a row [probability vector].
%               'pdfs' - component distributions. A cell, each element is a
%                   pdf structure.
%       'tranform' - transformation of the input data before density
%           evaluation. This should be a [general function].
%
%   See also ML_ESTPDF ML_RND

%   02-Aug-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

if isfield(param,'transform')
    x = ml_evalfun(x,param.transform);
end
    
switch param.name
    case 'exp'
        y = expcdf(x,param.beta);
    case 'norm'
        y = normcdf(x,param.mu,param.sigma);
    case 'poiss'
        y = poisscdf(round(x),param.lamda);
    case 'gamma'
        y = gamcdf(x,param.alpha,param.beta);
    case 'mix'
        k = length(param.ps);
        for i=1:k
            densities(:,i) = ml_cdf(x,param.pdfs{i});
        end

        y = densities*param.ps';
    otherwise
        error('Unrecognized pdf name');
end

