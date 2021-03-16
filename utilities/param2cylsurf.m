function nuclear_shape_model = param2cylsurf( cellparam, options )
% Train nuclear shape model using 3D HeLa cell nuclei images

% Author: Tao Peng
%
% Copyright (C) 2011-2016 Murphy Lab
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

% July 23, 2012 R.F. Murphy Add debug code to show average shape
%
% July 26, 2012 I. Cao-Berg Added a statement where it returns an empty model
%               when the spline features are empty
%
% August 2, 2012 I. Cao-Berg Added debugging statement that prints a message
%               if spline features calculation is not successful
%
% Jan 30, 2013 I. Cao-Berg Display of plots will happen if debug and
%               display flags are both true
%
% Jan 1, 2013 I. Cao-Berg Updated method to use display according to
%               verbose flag
%
% Jan 6, 2013 I. Cao-Berg Fixed bug in verbose, debug and display flags
%
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
% Sep 15, 2015 X. Ruan delete repeated code.
%
% March 4, 2016 I. Cao-Berg Included check that continues computations if
% cellcodes are empty
%
% March 4, 2016 I. Cao-Berg Included check that continues computations if
% some cellcodes are missing

if ~exist('options', 'var')
    options = [];
end
%icaoberg 06/02/2013

options = ml_initparam(options, struct('verbose', true, ...
    'debug', false, ...
    'display', false));

allcoefs = [];
allknotsx = [];
allknotsy = [];
allheight = [];

%D. Sullivan 6/5/13
%Load the preprocessed nucimage features

if ischar(cellparam{1})
    for i = 1:length(cellparam)
        cellparam{i} = load(cellparam{i}, 'nuc');
    end
end

nuc_model_file = [savepath filesep 'nuc_model.mat'];
for i = 1:length(cellparam)
    disp(['Concatenating spline features from image ' num2str(i) ] );
    
    if isempty(cellparam{i})
        if options.verbose
            disp( 'Spline features missing. Continuing to next image.' );
        end
        continue;
    end
    
    try
        u = (cellparam{i}.nuc.spfeat.coefs(:,1) + cellparam{i}.nuc.spfeat.coefs(:,end))/2;
        cellparam{i}.nuc.spfeat.coefs(:,1) = u;
        cellparam{i}.nuc.spfeat.coefs(:,end) = [];
        allcoefs = [allcoefs cellparam{i}.nuc.spfeat.coefs(:)];
        allheight = [allheight;cellparam{i}.nuc.spfeat.height];
        % --------------------------------------- %
        % To verify the constant knots assumption %
        allknotsx = [allknotsx;cellparam{i}.nuc.spfeat.knots{1}];
        allknotsy = [allknotsy;cellparam{i}.nuc.spfeat.knots{2}];
        % --------------------------------------- %
    catch err
        warning('Unable to load parameterization or parameterization is missing');
        getReport( err )
    end
    
end
allcoefs = allcoefs';

disp( 'Train nuclear shape model...' ); tic;
f_height = ml_estpdf(allheight,struct('name','norm'));
f_coef = ml_estpdf(allcoefs,struct('name','mvn'));

delta = 2*pi/360.;

averagecoef = reshape(f_coef.mu,size(cellparam{1}.nuc.spfeat.coefs));
averagecoef(:,end+1) = averagecoef(:,1);
avginstance = cellparam{i}.nuc.spfeat;
avginstance.coefs = averagecoef;
H = round(f_height.mu);

Phi = -pi:delta:pi;
Z = 0:(1/H):1;
[Phi_grid, Z_grid] = meshgrid(Phi,Z);
mesh_data = [Z_grid(:), Phi_grid(:)]';
nucsurf = reshape(fnval(avginstance,mesh_data),[length(Z),length(Phi)]);

if options.display
    try
        figure
        plotcylsurf(nucsurf,delta);
    catch
        if verbose; disp('Unable to display figure'); end
    end
end

disp( 'Creating structure...' );
nuclear_shape_model = struct('name','spsurf',...
    'surface',struct('form','B-',...
    'nknots_phi',5,'constknot_phi',[.25 .375 .5 .625 .75],...
    'nknots_h',1,'constknot_h',0.5,...
    'number',[4 9],...
    'order',[3 4],...
    'stat',f_coef),...
    'height',struct('stat',f_height));

nuclear_shape_model.type = 'cylindrical_surface';