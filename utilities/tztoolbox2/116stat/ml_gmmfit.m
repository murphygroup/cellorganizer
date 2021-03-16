function gmmpdf = ml_gmmfit(x,param)
%ML_GMMFIT Fit Gaussin mixture.
%   GMMPDF = ML_GMMFIT(X) returns a [pdf] representing a Gaussian mixture
%   model estimated the [feature matrix] X.
%   
%   GMMPDF = ML_GMMFIT(X,PARAM) specifies how to estimate the Gaussian
%   mixture by PARAM, which is a structure with the following fields:
%       'method' -
%           'mml' : a method based on minimal message length criteion.
%               reference: "Unsupervised Learning of Finite Mixture
%               Models", PAMI 24(3), 2002. p.381
%           'em' : EM algorithm
%           'fix': already know the labels
%               'labels' - a [label vector]. labels of X.
%   
%   See also

%   15-Aug-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('method','mml'));
gmmpdf.name = 'mix';

pdfname = ml_getgausspdfname(size(x,2));

switch param.method
    case 'mml'
        param = ml_initparam(param,struct('kmin',1,'kmax',10, ...
            'regularize',0,'th',1e-5,'covoption',0));
        [bestk,bestpp,bestmu,bestcov,dl,countf] = ...
            ml_mixtures4(x',param.kmin,param.kmax,param.regularize, ...
            param.th,param.covoption);
        
        if size(x,2)==1
            bestcov = bestcov.^0.5;
        end
        
        for i=1:bestk            
            gmmpdf.pdfs{i} = struct('name',pdfname, ...
                'mu',bestmu(:,i)','sigma',bestcov(:,:,i));          
        end
        gmmpdf.ps = bestpp;
    case 'em'
        param = ml_initparam(param,struct( ...
            'ncomp',3,'covartype','full','options',zeros(1,14),'weights',[]));
        mix = gmm(size(x,2),param.ncomp,param.covartype);    
        mix = gmminit(mix,x,param.options);
        mix = ml_gmmem(mix, x, param.options, param.weights);
        pdfname = ml_getgausspdfname(mix.ncentres);
        if size(x,2)==1
            mix.covars = mix.covars.^0.5;
        end
        for i=1:mix.ncentres
            gmmpdf.pdfs{i} = struct('name',pdfname, ...
                'mu',mix.centres(i,:),'sigma',mix.covars(:,:,i));
        end
        gmmpdf.ps = mix.priors;
    case 'fix' %already know the labels
        param = ml_initparam(param,struct('labels',[],'ispooled',0));
        if isempty(param.labels)
            error('You must sepcify labels to learn mixture models by labels');
        end
        
        if ~param.ispooled
            f.name = ml_getgausspdfname(size(x,2));
            nclass = max(param.labels);
            ns = [];
            gmmpdf.pdfs = {};
            for i=1:nclass
                y = x(param.labels==i,:);
%             disp(['data size: ' num2str(size(y,1))]);               
                if ~isempty(y)
                    if ~isempty(y)
                        constidx = ml_constidx(y);
                    end
                    if ~isempty(constidx)
                        y = ml_rmcol(y,constidx);
                    end
                    if isempty(y)
                        gmmpdf.pdfs{end+1} = struct([]);
                        ns(end+1) = 0;
                    else
                        ns(end+1) = size(y,1);
                        gmmpdf.pdfs{end+1} = ml_estpdf(y,f,param);
                        if ~isempty(constidx)
                            gmmpdf.pdfs{end}.transform = struct( ...
                                'funname','ml_rmcol','param',constidx);
                        end
                    end
                else
                    warning('a group is empty');
                end
            end
            gmmpdf.ps = ns/sum(ns);
        else
            [fsigma,mu,ps] = ml_estpoolcov(x,param.labels);
            gmmpdf = struct('name','mix');
            ncluster = max(param.labels);
            emptyidx = [];
            for i=1:ncluster
                idx = find(param.labels==i);
                if isempty(idx)
                    emptyidx = [emptyidx i];
                end
                gmmpdf.pdfs{i} = ...
                    struct('name','mvn','mu',mu(i,:),'sigma',fsigma);
            end
            gmmpdf.ps = ps;
            gmmpdf.pdfs(emptyidx) = [];
            gmmpdf.ps(emptyidx) = [];
        end
        
    otherwise
        error('Unrecognized method');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name=ml_getgausspdfname(n)

if n==1
    name = 'norm';
else
    name = 'mvn';
end

    
