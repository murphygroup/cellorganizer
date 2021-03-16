function model2info( model, fileID )
% MODEL2INFO List basic information about an SLML model.
%
% model   model structure (required)
% options (optional) options for modifying information displayed

% Ivan E. Cao-Berg
%
% Copyright (C) 2018 Murphy Lab
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

% 2018/03/18 icaoberg Updated method to save results to a file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'name' )
    fprintf( fileID,'%%%%%% Model name\n' );
    fprintf( fileID, ['''' model.name '''' ';\n'] );
else
    fprintf( fileID,'%%%%%% Model name\n' );
    fprintf( fileID, ['''UNSET''' '\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'id' )
    fprintf( fileID,'%%%%%% Model ID\n' );
    fprintf( fileID, ['''' model.id '''' ';\n'] );
else
    fprintf( fileID,'%%%%%% Model ID\n' );
    fprintf( fileID, ['''UNSET''' '\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIMENSIONALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'dimensionality' )
    fprintf( fileID,'%%%%%% Dimensionality\n' );
    fprintf( fileID, ['''' model.dimensionality ''''  ';\n'] );
else
    fprintf( 'Dimensionality: UNSET\n' );
    warning( 'This is an invalid model' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOCUMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'documentation' )
    fprintf( fileID,'%%%%%% Documentation\n' );
    fnames = fieldnames( model.documentation );
    fprintf( fileID,'%% Fields and values present in the documentation structure\n\n' );
    for i=1:1:length(fnames)
        fname = fnames{i};
        fprintf( fileID, ['%% ' fname '\n'] );
        fprintf( fileID, ['''' eval([ 'model.documentation.' fname]) ''''  ';\n'] );
        clear fname
    end
    clear fnames
end

try
    x=model.nuclearShapeModel;
    fprintf( fileID,'\n%%%%%% Nuclear shape model information\n' );
    model2info2(model.nuclearShapeModel, fileID )
catch
    fprintf( fileID,'\n%%%%%% Nuclear shape model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.cellShapeModel;
    disp(' ');
    fprintf( fileID,'\n%%%%%% Cell shape model information\n' );
    model2info2(model.cellShapeModel, fileID )
catch
    fprintf( fileID,'\n%%%%%% Cell shape model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.proteinModel;
    disp(' ');
    fprintf( fileID,'\n%%%%%% Protein pattern model information\n' );
    model2info2(model.proteinModel, fileID );
catch
    fprintf( fileID,'\n%%%%%% Protein pattern model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS DIFFEOMORPHIC? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_diffeomorphic( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    f = figure('visible','off');
    nimgs = size(model.cellShapeModel.positions,1);
    nlabels = size(model.cellShapeModel.positions,1);
    labels = reshape(repmat([1:nlabels],ceil(nimgs/nlabels),1),[],1);
    options.plot_dims = [1,2];
    options.subsize = 1000;
    showShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    
    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_shape_space.png" ><img src="show_shape_space_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 1. Shape space for input model. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
    
    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS PCA FRAMEWORK? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_pca_framework( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    f = figure('visible','off');
    showPCAShapeSpaceFigure(model);
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    
    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_shape_space.png" ><img src="show_shape_space_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 1. Shape space for input model. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
    
    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS CLASSIC MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_classic_model( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    %f = figure('visible','off');
    fprintf( fileID,'%% There are no additional figures for this model at the moment.\n' );
end
end

function model2info2( model, fileID )

fprintf( fileID,'%% Fields and values present in the submodel.\n\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'name' )
    fprintf( fileID, ['\n%% model.name\n'] );
    fprintf( fileID, ['''' model.name '''' ';\n'] );
else
    fprintf( fileID, ['%% model.name\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'id' )
    fprintf( fileID, ['%% model.id\n'] );
    fprintf( fileID, ['''' model.id '''' ';\n'] );
else
    fprintf( fileID, ['%% model.id\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'class' )
    fprintf( fileID, ['%% model.class\n'] );
    fprintf( fileID, ['''' model.class '''' ';\n'] );
else
    fprintf( fileID, ['%% model.class\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'type' )
    fprintf( fileID, ['%% model.type\n'] );
    fprintf( fileID, ['''' model.type '''' ';\n'] );
else
    fprintf( fileID, ['%% model.type\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'resolution' )
    fprintf( fileID, ['%% model.resolution\n'] );
    fprintf( fileID, ['''' mat2str(model.resolution) '''' ';\n'] );
else
    fprintf( fileID, ['%% model.resolution\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end
end
