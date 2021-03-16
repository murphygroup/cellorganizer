function f2 = ml_estpdf(x,f,param)
%ML_ESTPDF Estimate the distribution from data.
%   F2 = ML_ESTPDF(X,F) returns a [pdf], which is the estimation of the
%   [partial pdf] F on data X. If it is a univariate distribution, X must
%   be a column vector. If it is a multivariate distribution, X must be a
%   [feature matrix].
%   Notice: Although F.transform is usually a [general function], there is
%   one exception. F.transform will do PCA transformation when the name is
%   '_pca'. The number of PCA components can be specified by the field
%   'ncomp' in F.transform.param.
%   
%   F2 = ML_ESTPDF(X,F,PARAM) specifies how to estimate the distribution.
%   PARAM is a structure. Currently it is only available for the 'mvn' pdf.
%   If F.name is 'mvn', then PARAM has the following field:
%       'weights' - weights of the samples. It should have same number of 
%            rows as that of X.
%       'tz_estcov' - the parameter for ML_ESTCOV, which is used to
%           estimate the covariance matrix. This is only useful when the
%           [pdf] F has unknown covariance matrix. But the 'weights'
%           field in 'tz_estcov' has no effect.
%
%   See also ML_PDF ML_RND

%   03-Aug-2006 Initial write T. Zhao
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
    error('At least 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

ml_estpdf_checkvar(size(x,2),f.name);

if isfield(f,'transform')
    if strcmp(f.transform.funname,'_pca')
        f.transform.funname = 'ml_pcatrans';
        f.transform.param.basevec = princomp(x);
        if isfield(f.transform,'param')
            if isfield(f.transform.param,'ncomp')
                f.transform.param.basevec = ...
                    f.transform.param.basevec(:,1:f.transform.param.ncomp);
                rmfield(f.transform.param,'ncomp');
            end
        end
        f.transform.param.offset = mean(x,1);
        f.mu = zeros(1,size(f.transform.param.basevec,2));
    end
    x = ml_evalfun(x,f.transform);
end

switch f.name
    case 'trunc' %truncated function
        f = ml_initparam(f,struct('a1',min(x),'a2',max(x)));
        param = ml_initparam(param,struct('iter',100,'isshow',0));
        if f.a1>f.a2
            error(['The left bound should ' ...
                   'not be greater than the right round']);
        end
                   
        y = x;
        f2 = f;
        fs = {};
        lks = [];
        for k=1:param.iter
            comppdf = ml_estpdf(y,f.comppdf);
            f2.comppdf = comppdf; %intermediate results
            fs{k} = f2;
            
            cdf1 = ml_cdf(f.a2,comppdf) - ml_cdf(f.a1,comppdf);
            n = round(length(x)/cdf1);
            y = ml_rnd(comppdf,n);
            
            y(y>f.a1 & y<f.a2) = [];
            if isempty(y) %no new data points
                break;
            end
            y = [y;x];
            lks = [lks ml_loglk(x,f2)];
            
            if param.isshow==1
                subplot(1,2,1)
                ml_histpdfplot(x,f2);
                subplot(1,2,2)           
                plot(lks)
                drawnow
            end
        end
        if ~isempty(lks)
            [maxlk,idx] = max(lks);
            f = fs{idx(1)};
        else
            f = fs{1};
        end
    case 'tnorm' %truncated normal distribution
        y = x;
        a1 = min(x);
        for k=1:50
            f = ml_estpdf(y,struct('name','norm')); 
            cdf1 = normcdf(a1,f.mu,f.sigma);
            n = round(length(x)/(1-cdf1));
            y = ml_rnd(f,struct('n',n));
            y(y>a1) = [];
            y = [y;x];
            ml_histpdfplot(y,f);
            drawnow;
        end     
        f.name = 'tnorm';
        f.a1 = a1;
    case 'norm' %normal distribution       
        if ~isfield(f,'mu') & ~isfield(f,'sigma')
            [f.mu,f.sigma] = normfit(x);
        else
            if ~isfield(f,'mu')
                f.mu = mean(x);
            end

            if ~isfield(f,'sigma')
                f.sigma = sqrt(sum((x-f.mu).^2)/(length(x)-1));
            end
        end
    case 'mvn' %multivariate norml distribution
        param = ml_initparam(param, struct('tz_estcov', ...
                struct('method','mle'),'weights',[]));
       
        if isfield(param.tz_estcov,'weights')
            warning(['The field ''weight'' in param.tz_estcov has no ' ...
                     'effect']);
            param.tz_estcov.weights = param.weights;
        end
        
        if ~isfield(f,'mu')
            if isempty(param.weights)
                f.mu = mean(x,1);
            else
                f.mu = weights'*x/sum(weights);
            end
            
        end
        
        if ~isfield(f,'mu') & ~isfield(f,'sigma')
            f.sigma = ml_estcov(x,param.tz_estcov);
        else
            if ~isfield(f,'sigma')
                param.tz_estcov.mu = f.mu;
                f.sigma = ml_estcov(x,param.tz_estcov);
                %f.sigma = ml_cov(x,0,[],f.mu);
            end
        end    
    case 'gamma'
        parmhat = gamfit(x);
        f.alpha = parmhat(1);
        f.beta = parmhat(2);
    case 'exp'
        f.beta = expfit(x);
    case 'hist'
        [f.hist,occ1,occ2] = unique(x,'rows');
        f.hist(:,end+1) = ml_countnum(occ2)'/size(x,1);
    case 'mix'
        switch param.mixname
            case 'gmm' %gaussian mixture
                f = ml_gmmfit(x,param);
            otherwise
                error('Unrecognized mixture model name');
        end
    case 'bicond'
        f.pdf1 = ml_estpdf(x(:,1),f.pdf1);
        f.pdf2 = ml_estpdf(ml_evalfun({x(:,2),x(:,1)},f.relation),f.pdf2);
    otherwise
        error(['Unrecognized pdf:' f.name]);
end

f2 = f;


function ml_estpdf_checkvar(n,name)

if ~isempty(strmatch(name,{'norm','exp','gamma'}))
    if n>1
        error(['The pdf ' name ' does not support multivariate distributions.']);
    end
end

