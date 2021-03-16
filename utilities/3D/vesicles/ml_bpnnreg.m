function regmodel = ml_bpnnreg(x,y,t)

%ML_BPNNREG Multivariate regression by BPNN
%   REGMODEL=ML_BPNNREG(X,Y,T) returns a structure restoring 
%   the regression model of f(X)=E(Y|X). T is the structure of
%   parameters. It has the following fields:
%       norm, stop, randtrainsel: see ML_REGRESS
%       hidden, epochs, bp: see ML_MLPTRAINTEST
%
%   SEE ALSO ML_EVALBPNNREG

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%   16-May-2005 Initial write T. Zhao

if ~exist('t','var')
    t.norm=1;
end

if ~isfield(t,'norm')
    t.norm=1;
end

if ~isfield(t,'hidden')
    t.hidden = 20;
end

if ~isfield(t,'epochs')
    t.epochs=300;
end

if ~isfield(t,'bp')
    t.bp=1000;
end

if ~isfield(t,'stop')
    t.stop=1;
end

if t.stop==1
    t.norm=2;
end

trainset=[];
stopset=[];
rlabel=[];
slabel=[];
prep.zmean=[];
prep.zsdev=[];
prep.zscore=0;
constidx=ml_constidx(x);
if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end

%normalize upon all training data
if t.norm==1
    [x,zmean,zsdev]=ml_zscore(x);
    prep.zscore=1;
    prep.zmean=zmean;
    prep.zsdev=zsdev;
end

if t.stop==1
    nsample=size(x,1);
    
    if ~isfield(t,'randtrainsel')
        t.randtrainsel=[];
    end
    if isempty(t.randtrainsel)
        t.randtrainsel=randperm(nsample);
    end
    
    stopsize=round(nsample/3);
    stopset=x(t.randtrainsel(1:stopsize),:);
    trainset=x(t.randtrainsel(stopsize+1:end),:);
    tmpytr=y(t.randtrainsel(stopsize+1:end),:);
    tmpyst=y(t.randtrainsel(1:stopsize),:);
    
    tmptrainset=trainset;
    
    switch(t.norm)
    case 1
        [trainset,zmean,zsdev]=ml_zscore(trainset);
        stopset=ml_zscore(stopset,zmean,zsdev);
        prep.zscore=1;
        prep.zmean=zmean;
        prep.zsdev=zsdev;
    case 2
        [trainset,zmean,zsdev]=ml_zscore(tmptrainset);
        stopset=ml_zscore(stopset,zmean,zsdev);
        prep.zscore=1;
        prep.zmean=zmean;
        prep.zsdev=zsdev;
    end
    
    % Train the network using the training and stop data
    [trainnetout stopnetout imin net] = ...
        ml_mlptraintest(trainset, tmpytr,stopset, tmpyst, t.hidden, t.epochs,t.bp) ;
else
    net = mlp(size(x,2), t.hidden, size(y,2), 'logistic') ;
    
    roptions = zeros(1,18) ;
    roptions(1) = 1 ;   % Output sse values
    %roptions(1) = -1 ;  % Output nothing 
    roptions(14) = t.epochs ;  % Number of epochs (train one epoch at a time)
    roptions(17) = 0.9 ;
    roptions(18) = 0.001 ;
       
    [net, roptions] = netopt(net, roptions, x, y, 'graddesc') ;     
end

regmodel=struct('modelname','bpnn','modeltype',...
    'network','trained',net,'t',t,'prep',prep);

%
% Summarize network performance
% y = mlpfwd(net, sample);
