function y = ml_pdf(x,param)
%ML_PDF Computes a probablity density or mass function.
%   Y = ML_PDF(X,PARAM) returns the probablity or probablity density of the
%   data X. If the pdf is univariate, then each element of X is a data
%   point. If the pdf is multivariate, then each row is a data point. But
%   there are exceptions. For finite mixture distribution, each row of X
%   is a data point no matter whether it is univaraite or multivaraite.
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
%           'mvn' : multivariate normal distribution
%               'mu' -mean, a row vector
%               'sigma' - covariance matrix 
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
%           'bicond' : bivariate conditional distribution. X must have 2 
%                columns. The probability function is 
%                P(X)=P(X(:,1))*P(f(X)). Therefore PARAM has the
%                following fields:
%                    pdf1 - [pdf] for X(:,1)
%                    relation - [general function] of {X(:,2),X(:,1)}. It must 
%                        be a bivariate function.
%                    pdf2 - [pdf] for f(X(:,2),X(:,1))
%           'trunc' : truncated distribution
%                comppdf - pdf becfore truncation
%                a1 - left bound
%                a2 - right bound
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
        y = exppdf(x,param.beta);
    case 'norm'
        y = normpdf(x,param.mu,param.sigma);
    case 'tnorm'
        y = zeros(length(x),1);
        cdf1 = normcdf(param.a1,param.mu,param.sigma);
        y(x>=param.a1) = normpdf(x(x>=param.a1),param.mu,param.sigma)/(1-cdf1);
    case 'mvn'
        if all(param.sigma(:)==0)
            mus = repmat(param.mu,size(x,1),1);
            y = double(sum(~(x==mus),2)==0);
        else
            y = mvnpdf(x,param.mu,param.sigma);
        end        
    case 'poiss'
        y = poisspdf(round(x),param.lamda);
    case 'gamma'
        y = gampdf(x,param.alpha,param.beta);
    case 'mix'
        k = length(param.ps);
        for i=1:k
	     if isempty(param.pdfs{i})                %%jl++, 10-31-2011
		 densities(:,i) = 0;
	     else
            %densities(:,i) = ml_pdf(x,param.pdfs{i});  %%jl--, 10-31-2011
            densities(:,i) = prod(ml_pdf(x,param.pdfs{i}));  %%jl++, 10-31-2011
	     end
        end

        y = densities*param.ps';
    case 'bicond'
        y = ml_pdf(x(:,1),param.pdf1);
        y = y.*ml_pdf(ml_evalfun({x(:,2),x(:,1)},param.relation), ...
                      param.pdf2);
    case 'unnorm'
        y = ml_evalfun(x,param.f);
    case 'trunc'
        y = zeros(size(x));
        cdf1 = ml_cdf(param.a2,param.comppdf) - ...
            ml_cdf(param.a1,param.comppdf);
        idx  = find(x>=param.a1 & x<=param.a2); 
        y(idx) = ml_pdf(x(idx),param.comppdf)/cdf1;
    otherwise
        error('Unrecognized pdf name');
end

