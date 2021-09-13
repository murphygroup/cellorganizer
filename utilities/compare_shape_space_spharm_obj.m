
function [ consolidated_model ] = compare_shape_space_spharm_obj( models,fileID,options )

 consolidated_model=consolidateshapespacemodels(models{1},models{2});

 f = figure('visible','off');
 options.viewangle = [0,90]; %down z axis
 updateLabels(consolidated_model,options);
 saveas( f, sprintf('show_shape_space.png', 'png'));
 I = imread( 'show_shape_space.png' );
 I = imresize( I, 0.50 );
 imwrite( I, 'show_shape_space_thumbnail.png' );
 img2html(fileID,'show_shape_space.png','show_shape_space_thumbnail.png','Joint shape space from SPRM-RPDM models.');
end

