function [lk,lks] = ml_loglk(x,f,param)
%ML_LOGLK Calculate log likelihood.
%   LK = ML_LOGLK(X,F) calculate log likelihood of data X with the [pdf] F.
%   
%   LK = ML_LOGLK(X,F,WEIGHTS) assumes the data are weighted. WEIGHTS is a
%   column vector that has the same length as the number of rows in X. So
%   the ith element of WEIGHTS is the weight of the ith row of X.
%
%   [LK,LKS] = ML_LOGLK(X,F,WEIGHTS) also returns likelihood for each data
%   point.
%
%   See also

%   24-Aug-2006 Initial write T. Zhao
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

if ~exist('weights','var')
    weights = [];
end

if isfield(f,'transform')
    x = ml_evalfun(x,f.transform);
    rmfield(f,'transform');
end

switch f.name
    case 'norm'
        lk = 0.5*(x-f.mu).^2/f.sigma^2; 
        if isempty(weights)
            total = length(x);
        else
            total = sum(weights);
            lk = lk.*weights;
        end   
        lk = -0.5*log(2*pi*f.sigma.^2) - lk;
    case 'mvn'
        x = ml_addrow(x,-f.mu);
        r = chol(f.sigma);
        x = x/r;
        lk = 0.5*sum(x.^2,2);
        
        if isempty(weights)
            total = size(x,1);
        else
            total = sum(weights);
            lk = lk.*weights;
        end  
        
        lk = -0.5*size(x,2)*log(2*pi)-0.5*log(det(f.sigma))-lk;
    otherwise
        lk = log(ml_pdf(x,f));    
        %Estimate lower bound if the pdf is out of floating range
        if any(isinf(lk)) & strcmp(f.name,'mix')
            idx = find(isinf(lk));
            x = x(idx,:);
            lk(idx) = 0;
            tmplks = [];
            for i=1:length(f.pdfs)
                [tmplk,tmplks(:,i)] = ml_loglk(x,f.pdfs{i});
                tmplks(:,i) = tmplks(:,i)+log(f.ps(i));
                lk(idx) = max(tmplks,[],2);
            end
%             if any(isinf(lk(idx)))
%                 keyboard;
%             end
            
        end
        if ~isempty(weights)
            lk = lk.*weights;
        end
end

if nargout==2
    lks = lk;
end

lk = sum(lk);
