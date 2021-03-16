function proteinmodel = param2protein_gaussian_model( cellparam, options)
% Train objectized protein and position model using 3D HeLa GFP images
%
% Adapted from train_protein_model2 - grj 10/13/15

%
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

options = ml_initparam(options, struct('verbose', false, 'debug', false));

for i = 1:length(cellparam)
    if ischar(cellparam{i})
        cellparam{i} = load(cellparam{i}, 'prot');
    end
end

if iscell(cellparam(1))
    cellparam = [cellparam{:}];
end

protparam = [cellparam.prot];

% 
% disp( 'Compiling Gaussian mixtures...' );
% 
% tic; 
% 
% 
% compile_gaussobjs( options.savefitdir, options.tempparent, options )
% 
% toc

disp( 'Learning object models' )
tic
proteinmodel.size = gmm_objsizefit(vertcat(protparam.gauss_objsize));
toc

disp( 'Learning frequency model' )
tic
%D. Sullivan 6/12/13 - changed argument to param struct.
% model.frequency = getProteinFrequencyModel( savedir );
proteinmodel.frequency = getProteinFrequencyModel( cellfun(@length, {protparam.gauss_intens}) );
toc

disp( 'Learning object position models' )
tic
pos = [protparam.pos];
proteinmodel.position.beta = vertcat(pos.beta);
toc

%D. Sullivan 2/22/13 also need to save the resolutions of the models 
%D. Sullivan 6/12/13 added adjustment for resolution formerly dealt with in
%gmm_objempdistr.m. Now setting this in per-cell calculations is not
%possible due to parallelization, but computation is still done.

if isfield(options.model,'proteinUpsampleZ')
    options.model.protein_resolution(3) = options.model.protein_resolution(3)/options.model.proteinUpsampleZ;
end
proteinmodel.resolution = options.model.protein_resolution;
proteinmodel.class = 'object';

end

function model = gmm_objsizefit( objsize )
% Learn object size model for each location pattern

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 16, 2012 I. Cao-Berg Removed protein type from function call
% July 30, 2012 I. Cao-Berg Updated object size model to be a fit 
%                           from an exponential to a lognormal
% August 6, 2012 D. Sullivan major bug fix: added separate array for xbin 
%           values to prevent the x values from being overwritten before 
%           the z_x distribution is learned
% August 7, 2012 R.F. Murphy fix loop limit
%  
%
% Edited: D. Sullivan 8/6/2012 undid some of the previous changes such as
% binsizes. Go from 0-xobjsize==30 (arbitrarily chosen) and use a max_xbin
% = 15 (arbitrarily chosen).
%
% Edited: D. Sullivan 8/7/2012 reverted entirely to original version to fix
% issue with objects no longer appearing gaussian. Then increased range of
% x from 0-xobjsize==30 (arbitrarily chosen) and use a max_xbin
% = 15 (arbitrarily chosen).
%
%D. Sullivan 6/11/13 - got rid of all the resolution dependencies in the
%                      size restriction. based on 95% confidence interval
%                      now. Also fitting x only based on 95% confidence of
%                      object sizes now instead of whole distribution
%GRJ 10/13/15 - moved to param2protein_gaussian_model.m


%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval
%arbitrarily chosen to be the maximum size of xsigma to consider
% max_xbin = 15;

%D. Sullivan 6/12/13 added loop to gather per-cell computed objects 
% files = ml_dir([datadir filesep '*.mat']);
% objsizetot = [];
% tic
% for i = 1:length(files)
% load([datadir filesep files{i}]);
% objsizetot = [objsizetot;real(objsize)];
% end
% objsize = objsizetot;
% clear objsizetot
% toc

x = objsize(:,1);
y = objsize(:,2);
z = objsize(:,3);

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval X DIMENSION
tmpx = sort(x);
x95 = tmpx(1:floor(0.95*length(x)));
max_xbin = ceil(x95(end));

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval Y DIMENSION
ty = sort(y);
y95 = ty(1:floor(0.95*length(ty)));
optvars.alpha1y = y95(end);%(LAMP this = ~5.9)vs6.5
optvars.alpha2y = std(y95);%(LAMP this = ~1.3)vs1.4

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval Z DIMENSION
tz = sort(z);
z95 = tz(1:floor(0.95*length(tz)));
optvars.alpha1z = z95(end);
optvars.alpha2z = std(z95);

%icaoberg 7/30/2012
%lambda = expfit(x);
%model.x = struct('name','exp','beta',lambda);
%D. Sullivan, 6/11/13 now only fitting 95% confidence sizes
[parmhat,parmci] = lognfit( x95 );
model.x = struct('name', 'lognorm', 'mu', parmhat(1), 'sigma', parmhat(2) );

% Fit the y-x conditional normal distribution
Y_X = cell(max_xbin,1);
for i=1:max_xbin
    y_given_x = y(x>i-1&x<=i);
    Y_X{i} = y_given_x;
    %D. Sullivan 6/11/13 - added else in case there are no observed objects
    %in that bin
    if ~isempty(y_given_x)
        [mu(i),sigma(i)] = normfit(y_given_x);
    else
        
    end
end

% dpsulliv 8/7/12 added max_xbin as a variable
% x = .5:14.5;
xbin = .5:max_xbin-0.5;

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval
% optvars.alpha1 = 6.5;
% optvars.alpha2 = 1.4;
% dpsulliv 8/7/12 added max_xbin as a variable
% u = mu(1:15);
u = mu(1:max_xbin);
xbin(u>=optvars.alpha1y) = [];
u(u>=optvars.alpha1y) = [];
optvars.beta1y = sum(-xbin.*log(1-u/optvars.alpha1y))/sum(xbin.^2);
% dpsulliv 8/7/12 added max_xbin as a variable
% u = sigma(1:15);
u = sigma(1:max_xbin);
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
xbin(u>=optvars.alpha2y) = [];
u(u>=optvars.alpha2y) = [];
optvars.beta2y = sum(-xbin.*log(1-u/optvars.alpha2y))/sum(xbin.^2);

model.y_x = struct('name','norm',...
                    'a1',optvars.alpha1y,...
                    'b1',optvars.beta1y,...
                    'a2',optvars.alpha2y,...
                    'b2',optvars.beta2y);

% Fit the z-x conditional normal distribution
%clear mu sigma
% dpsulliv 8/7/12 expanded the range to 30 
% Z_X = cell(20,1);
% for i=1:20
Z_X = cell(max_xbin,1);
for i=1:max_xbin
    z_given_x = z(x>i-1&x<=i);
    Z_X{i} = z_given_x;
    if ~isempty(z_given_x)
        [mu(i),sigma(i)] = normfit(z_given_x);
    end
end

% dpsulliv 8/7/12 added max_xbin as a variable
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval
% optvars.alpha1 = 3.5;
% optvars.alpha2 = 1.4;
% dpsulliv 8/7/12 added max_xbin as a variable
% u = mu(1:15);
u = mu(1:max_xbin);
xbin(u>=optvars.alpha1z) = [];
u(u>=optvars.alpha1z) = [];
optvars.beta1z = sum(-xbin.*log(1-u/optvars.alpha1z))/sum(xbin.^2);
% dpsulliv 8/7/12 added max_xbin as a variable
% u = sigma(1:15);
u = sigma(1:max_xbin);
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
xbin(u>=optvars.alpha2z) = [];
u(u>=optvars.alpha2z) = [];
optvars.beta2z = sum(-xbin.*log(1-u/optvars.alpha2z))/sum(xbin.^2);

model.z_x = struct('name','norm',...
                    'a1',optvars.alpha1z,...
                    'b1',optvars.beta1z,...
                    'a2',optvars.alpha2z,...
                    'b2',optvars.beta2z);
end

function frequency = getProteinFrequencyModel( numberOfObjects )
    %GETPROTEINFREQUENCYMODEL Returns the protein frequency model component
    %using the intermediate results calculated by CellOrganizer.

    % Author: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
    %
    % Edited: D. Sullivan 6/12/13 - refactored the code. now reading objects
    %                               in from the object_stats folder. 
    %                               (contains object intensities and sizes,
    %                               + totnumber for each cell)



    % end

    numberOfObjects=numberOfObjects(find(numberOfObjects~=0));

    [parmhat,parmci] = lognfit( numberOfObjects );
    frequency = [];
    frequency.mu = parmhat(1);
    frequency.sigma = parmhat(2);
end

function beta = objpos_run2( tempfilepath, param )
% Learns the protein object position model

% Tao Peng
%
% Copyright (C) 2012-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 17, 2012 I. Cao-Berg Added parameter structure to the method and fixed a
%                     bug where the method complained because of a nonexisting
%                     preallocated beta
% Feb 22, 2013 D. Sullivan Changed downsampling to be defined by vector
%                          rather than hard coded
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
% June 12, 2013 D. Sullivan refactored to run with per-cell precomputed
%                           positions
% July 25, 2013 G. Johnson  Remove NANs from object positions so that
%                           'beta' is always legal
%

beta = [];
if nargin == 4
    param = [];
    verbose = false;
    debug = false;
elseif nargin > 5
    error('Wrong number of input arguments');
else
    try
        verbose = param.verbose;
        if ~islogical( verbose )
            verbose = false;
        end
    catch
        verbose = false;
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

%D. Sullivan 6/12/13 don't need this anymore due to refactoring
%D. Sullivan 2/22/13
%changed from being hard coded downsample to param based downsampling
% downsample = param.downsampling;
% downsample = [5 5 1];

%icaoberg 15/5/2013
% dna_image_files = ml_ls( dnaImagesDirectoryPath );
% cell_image_files = ml_ls( cellImagesDirectoryPath );
% prot_image_files = ml_ls( protImagesDirectoryPath );
% try
%     masks_image_files = ml_ls( param.masks );
% catch
%     mask_image_files = '';
% end

if ~exist( [tempfilepath filesep 'beta_all.mat'] )
    Xtot = [];
    Ytot = [];
    files = ml_dir( [tempfilepath filesep '*_Beta.mat'] );
    imageID = [];
    for i = 1:1:length(files)
%         try
            disp(['Image ' num2str(i)])
            load([tempfilepath filesep files{i}]);
            
            %grj 7/23/13
            rminds = any(isnan(X),2) | any(isnan(Y),2);
            
            X(rminds,:) = [];
            Y(rminds,:) = [];
            
            Xtot = [Xtot;X];
            Ytot = [Ytot;Y];
%             betatot(size(betatot,1)+1,:) = ml_logreg(x,distcodes(:,3))';
%             imageID = [imageID;repmat(size(betatot,1)+1,size(x,1),1)];
%         catch
%             disp(['Ignoring image:' num2str(i)]);
%         end
    end
    
    beta = ml_logreg(Xtot,Ytot)';
    save([tempfilepath filesep 'beta_all.mat'],'beta','X','Y','imageID')
else
    load( [tempfilepath filesep 'beta_all.mat'] );
end

end