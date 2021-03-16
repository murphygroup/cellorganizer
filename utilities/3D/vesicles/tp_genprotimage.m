function [protimg] = tp_genprotimage( model, param )
% Generate protein channel using a protein model, giving synthesized
% nuclear and cell filled images.

% Author: Tao Peng
%
% Copyright (C) 2008-2013 Murphy Lab
% Carnegie Mellon University
%
% February 29, 2012 D. Sullivan
% March 8, 2012 I. Cao-Berg Fixed bug in ml_celldistcode input arguments
% March 13, 2012 Devin Sullivan Moved gaussobj addition into tp_gengaussobjs
% March 14, 2012 Devin Sullivan Added sampling density
% March 18, 2012 R.F.Murphy Force number of objects to reasonable range;
% Must provide filled cell and nuclear image so that these are not recalculated for each
% pattern (via param.nucleus and param.cell)
% March 19, 2012 R.F.Murphy Eliminate bwperim calls to speed up (assumes
% param.cell gives same result as bwperim(param.cell)
% August 1, 2012 R.F.Murphy Narrow range on generated number of objects
% August 4, 2012 D. Sullivan Passed param.nucleus and param.cell from parent call
% to save memory
% November 15, 2012 I. Cao-Berg Changed warning to regular display
% March 14, 2013 D. Sullivan 3/14/13 added param structure to tp_gengaussobjs call
% July 18, 2013 D. Sullivan added pos to the outputs of tp_gengaussobjimg.m
%                           to allow for object motion model. Also check if
%                           'param.randomwalk' then save objects and
%                           positions
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

protimg=[];
if nargin ~= 2
  warning( 'CellOrganizer: Wrong number of input arguments' )
end

try
   verbose = param.verbose;
catch
   verbose = true;
end

try
   fileID = param.fileID;
catch
   fileID = [];
end

%nucimg = param.nucleus; %don't need these copies
%cellimg = param.cell;
sampling = param.sampling.method;

%D. Sullivan 12/4/14 - created support for multiple numbers to be specified
%If NaN is passed, the number is ignored and sampled from the distribution
% try
    if isfield(param,'numberOfGaussianObjects')
        if length(param.numberOfGaussianObjects)>1
            N = param.numberOfGaussianObjects(param.currentmodelnum);
        else
            N = param.numberOfGaussianObjects;
        end
    else
        N = [];
    end
% catch
%   N = [];
% end

% 28-OCT-2011 R.F. Murphy Use lognormal distribution model if N not specified
% D. Sullivan 12/4/14 - also use lognormal if NaN is specified
if isempty( N ) || isnan(N)
   %Sample randomly from the lognormal distribution to determine the
   %number of gaussian objects; trimmed to reasonable probability interval
   %to avoid unrealistically small or very large numbers of objects
   N=0; PofN=0;
   while ( N < 1 || PofN < 0.1 || PofN > 0.9 )
       N = ceil(lognrnd(model.frequency.mu,model.frequency.sigma,1,1));
       PofN = logncdf(N,model.frequency.mu,model.frequency.sigma);
   end
end

if verbose
   disp( ['Setting frequency model for protein shape model instance.']);
   disp( ['Number of objects set to ' num2str(N) '.' ] );
end

if verbose
   disp( 'Setting position model for protein shape model instance.' );
end

beta = model.position.beta;

if size(beta,1) > 1
    beta = beta(randi(size(beta,1)),:);
end

%devins 8/4/2012
%nucimg = imfill( nucimg>0, 'holes' );
%cellimg = imfill( cellimg>0, 'holes' );
%cellEdge = bwperim( param.cell, 18 );
%positions = (param.cell>0) - param.cell;
protimg = zeros(size(param.cell));

%distcodes = distances from nucleus and edge 
%coords = xyz coordinates for each point inside the cell and outside nucleus
%angles = polar coordinate angles for each point

try
    allowedCompartments = model.cytonuclearflag;
catch
    warning( 'No restriction on object position specified, defaulting to nucleus and cytoplasm' );
    allowedCompartments = 'all';
end

%icaoberg
if strcmpi( allowedCompartments, 'cyto' )
   allowedCompartments = 'cyt';
end

if strcmpi( allowedCompartments, 'nuclear' )
   allowedCompartments = 'nuc';
end

[distcodes,coords,angles] = ...
    ml_celldistcode2( param, allowedCompartments);
   %ml_celldistcode2( param, positions, 1, allowedCompartments);
nullidx = find(distcodes(:,1)==0 & distcodes(:,2)==0);
distcodes(nullidx,:)=[];
coords(nullidx,:)=[];
angles.theta(nullidx,:)=[];
angles.phi(nullidx,:)=[];

%normdists = fractional distance between nucleus and edge
normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);
x = ml_mappos([normdists angles.theta angles.phi]);

try
   e = x*beta;
catch
   beta = beta';
   e = x*beta;
end

%P = probability density associated at each position
P = exp(e)./(1+exp(e));

inds = sub2ind(size(param.cell), coords(:,1), coords(:,2), coords(:,3));

probimg = zeros(size(param.cell));
probimg(inds) = P;


if verbose
   disp( 'Sampling positions for all objects.' );
end


try
   objectmethod = param.sampling.method;
catch
   objectmethod = 'disc';
end

%what tells you the number of samples?
%set sampling density to 0 for cases where it will not be used
%(disc,trimmed)
samplingDensity = 0;
%if we are sampling, figure out the sampling density
if strcmpi(param.sampling.method,'sampled');
  try
    samplingDensity = param.sampling.density;
  catch 
    %icaoberg 11/15/2012
    %warning(['CellOrganizer: No sampling density specified with ',...
    %'sampling method set to sampled. Defaulting to 100']);
    disp(['No sampling density specified with ', ...
    'sampling method set to sampled. Defaulting to 100']);
    samplingDensity = 100; 
  end
end

disp('Generating objects');
%save('tempJustBeforeObjs.mat');

%D. Sullivan 3/14/13 added param structure to tp_gengaussobjs call
% [GaussObjects,protimgtmp] = tp_gengaussobjimg(model.size,N,...
%                 probimg,objectmethod,samplingDensity,param);
%D. Sullivan 7/18/13 added pos to outputs to track object positions for
%motion model
[GaussObjects,protimgtmp,pos] = tp_gengaussobjimg(model.size,N,...
    probimg,objectmethod,samplingDensity,param);

%D. Sullivan 7/18/13 saving Gaussian objects and their positions for later
%rendering if we are doing a walk.
% if isfield(param,'randomwalk') && param.randomwalk == true
%     save([param.targetdirectory filesep 'OriginalObjects.mat'],...
%         'GaussObjects','pos');
save([param.temporary_results filesep 'OriginalObjects.mat'],...
    'GaussObjects','pos');
% end
            
%eliminate samples outside the cell.
switch allowedCompartments
    case 'cyt'
        codemask = double(param.cell)-double(param.nucleus);
    case 'nuc'
        codemask = double(param.nucleus - bwperim(param.nucleus));
    case 'all'
        codemask = double(param.cell- bwperim(param.cell));
    otherwise
        error(['Unrecognized location name: ' loc]);
end

protimg = protimgtmp.*uint8(codemask);