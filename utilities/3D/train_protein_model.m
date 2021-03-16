function model = train_protein_model( ...
    dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, ...
    param )
% Train objectized protein and position model using 3D HeLa GFP images

% Tao Peng
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 17, 2012 I. Cao-Berg Added parameter structure to method
% April 24, 2012 R.F. Murphy Correction Gaussian object frequency model
% June 19, 2012 I. Cao-Berg Uncommented training of the frequency model
% Jan 1, 2013 I. Cao-Berg Updated method to use display according to
%               verbose flag
% Feb 22, 2013 D. Sullivan Added param structures to several training
%               functions so that models trained are adjusted and recorded
%               properly.
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
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

if nargin == 3
  param = [];
  verbose = true;
  debug = false;
else
  try 
   verbose = param.verbose;
   if ~islogical( verbose )
    verbose = true;
   end
  catch
   verbose = true;
  end

  try
   debug = param.debug;
   if ~islogical( debug )
    debug = false;
   end
  catch
    debug = false;
  end
end

%D. Sullivan 6/12/13 all temp dirs should be set up already in the param
%structure and used in the per-cell calculations.
% temporaryFilesFolder = [ pwd filesep 'temp' filesep 'protein_objects_gaussian' ];
% if ~exist( temporaryFilesFolder, 'dir' )
%     mkdir( temporaryFilesFolder );
% end
% 
% % Extract protein objects from 4 image sets (note this step can be parallelized)
% savedir = [ temporaryFilesFolder filesep 'original_objects'];
% if ~exist( savedir, 'dir' )
%     mkdir( savedir );
% end

%D. Sullivan 6/12/13 This is done per-cell now. 
% disp( 'Finding objects in images...' ); tic;
%D. Sullivan 2/22/13 added param input to findobjs_run to pass the
%resolution of the prot image. It also passes back a param struct with
%adjusted resolutions
% findobjs_run( imgdir, savedir);
% [param] = findobjs_run( ...
%     dnaImagesDirectoryPath, ...
%     cellImagesDirectoryPath, ...
%     proteinImagesDirectoryPath, savedir, param );
% toc

disp( 'Fitting Gaussian mixtures to the objects...' );
savedir = [ pwd filesep 'temp/protein_objects_gaussian/object_gaussians'];
if ~exist(savedir,'dir')
    mkdir(savedir);
end

tic; 
%D. Sullivan 2/22/13 added param input to learngaussobjs_run2 to pass the
%resolution of the prot image
% learngaussobjs_run2(imgdir,savedir);
%D. Sullivan 6/12/13 - already learned objects in a per-cell manner, now
%only need to compile them
compile_gaussobjs( param.savefitdir, param.tempparent, param )
% learngaussobjs_run2( ...
%     dnaImagesDirectoryPath, ...
%     cellImagesDirectoryPath, ...
%     proteinImagesDirectoryPath, savedir, param);
toc

if verbose; disp( 'Learning object models...' ); end;
gaussianobjdir = './temp/protein_objects_gaussian/object_gaussians';
savedir = './temp/protein_objects_gaussian/object_stats';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

if verbose; tic; end
%D. Sullivan 2/24/13 added param struct to track resolution and pass in the
%upsampling in the z dimension (param.proteinUpsampleZ)
param = gmm_objempdistr(gaussianobjdir,savedir,param);
%gmm_objempdistr(gaussianobjdir,savedir);
model.size = gmm_objsizefit(savedir);
if verbose; toc, end

if verbose; disp( 'Learning protein frequency model...' ); end
savedir = './temp/protein_objects_gaussian/object_gaussians';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

if verbose; tic; end
if verbose; disp( 'Learning frequency model...' ); end
model.frequency = getProteinFrequencyModel( savedir );
if verbose; toc, end

if verbose; disp( 'Learning object position models...' ); end
savedir = './temp/protein_objects_gaussian/object_positions';
if ~exist(savedir,'dir')
    mkdir(savedir)
end

tic;
%D. Sullivan 2/22/13 added param input to objpos_run to pass the
%resolution of the prot image
% model.position.beta = objpos_run(imgdir,savedir);
model.position.beta = objpos_run( ...
    dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, savedir, param );
toc

%D. Sullivan 2/22/13 also need to save the resolutions of the models 
model.resolution = param.model.protein_resolution;

% Save protein models
% disp( 'Saving protein models' );
% savedir = './inter_results/protein_objects_gaussian';
% save([savedir filesep 'model'],'model');
