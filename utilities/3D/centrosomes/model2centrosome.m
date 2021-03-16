function [centrosome,resolution] = model2centrosome(segdna, segcell, model,param)
%MODEL2CENTROSOME Generates a centrosome location from a generative model of protein location.

% Author: Jieyue Li
% Edited: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
% Edited: Devin Sullivan March, 18, 2012
%   - added up and down sampling 
%   - changed from model.parameters.slice to model.parameters.zloc, change
%   should be backwards compatible
% Edited: Devin Sullivan Feb, 27, 2013
%   Removed resolution adjustments, now done in model2instance.m
%
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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


% Dbothfin: cytosolic space, 0 - inside, 1 - outside;
% segdna: segmented nucleus, 0 - outside, 1 - inside;
% segcell: segmented cell,  0 - outside, 1 - inside.
% rand_seed: rand seed for generation.

centrosome = [];

 %If the model resolution is different than the cell resolution resize the
 %cell. - DPS 2015/04/14
 if any(param.resolution.cell~=model.resolution)
     segdna = AdjustResolutions(segdna,param.resolution.cell,model.resolution,1);
     segcell = AdjustResolutions(segcell,param.resolution.cell,model.resolution,1);
 end

Dbothfin = (logical(~segcell)+logical(segdna));

%D. Sullivan 2/27/13 Removed all resizing code, now done in
%model2instance.m

if ~isfield( model, 'class' )
   warning( 'CellOrganizer: Model class not present.' );
   return;
end

if ~strcmpi( lower(model.class), 'centrosome' )
   warning( 'CellOrganizer: Model class must be centrosome' );
   return;
end

if ~isfield( model, 'dimensionality' )
   warning( 'CellOrganizer: Model dimensionality not present. Assuming 3D.' );
end

if ~strcmpi( lower(model.dimensionality), '3d' )
  warning( ['CellOrganizer: Model dimensionality not 3D. This version only supports 3D ' ...
      'generative models of centrosomes.'] );
  return;
end

try
  cent_dist = model.parameters.distance;
catch
  warning( 'CellOrganizer: Incomplete model. Parameter distance not present or invalid.' );
  return;
end

try
  p = model.parameters.p;
catch
  warning( 'CellOrganizer: Incomplete model. Parameter p not present or invalid.' );
end

try
    try
      fractionZloc = model.parameters.zloc;
      slice = floor(fractionZloc*size(Dbothfin,3));
    catch
        warning('CellOrganizer: Old model type, using model.parameters.slice instead of model.parameters.zloc');
        slice = model.parameters.slice;
    end
  
catch
  warning( 'CellOrganizer: Incomplete model. Parameter slice not present or invalid.' );
end

try
 R = gamrnd(p(1),p(2));

 [D1] = bwdist(segdna(:,:,slice)>0);
 [D2] = bwdist((~segcell(:,:,slice))>0);
 
 Dbothfin = double(~Dbothfin(:,:,slice));
 Dbothfin(Dbothfin==0) = NaN;
 D1 = D1.*double(Dbothfin); 
 D2 = D2.*double(Dbothfin);
 
 DD = D1./D2;
 DD = DD - R;
 DD = abs(DD);%*4));%removed DPS 3/18/12 at the direction of Jieyue 
 

 [a,b] = find(DD==nanmin(nanmin(DD)));
 aaa = randperm(length(a));
 
 imgcent_coordinate = [a(aaa(1)),b(aaa(1)),slice];%slice*4];
 centrosome = imgcent_coordinate;
 
%D. Sullivan 2/27/13 Removed upsampling code, now done in model2instance.m
%also tracked resolution.
resolution = model.resolution;

%save temp
catch
 warning('CellOrganizer: Cannot generate centrosomal location');
end


