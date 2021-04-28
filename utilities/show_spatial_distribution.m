function [ output_args ] = show_spatial_distribution( model,fileID)
%Author : Serena Abraham



    f = figure('visible','off');
    spatial=model;
    [x1,y1,z1]=sph2cart(spatial.anglestheta,spatial.anglesphi,spatial.normdists);
    [X,Y,Z] = sphere;
    scatter3(X(:),Y(:),Z(:),'filled','MarkerFaceColor',[0.7 0.7 0.7]);
    hold on;
    scatter3(x1,y1,z1,'filled','b');
    legend('cell boundary','model');
    title('Spatial Distribution of objects');
    hold off;
    saveas( f, 'show_spatial_distribution.png', 'png');
    I = imread( 'show_spatial_distribution.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_spatial_distribution_thumbnail.png' );
    img2html(fileID,'show_spatial_distribution.png','show_spatial_distribution_thumbnail.png','Spatial Distribution across both models');
    
    
    
    %spatial positions fit
    %Generate all possible angles
    angle_range=linspace(0,6.28319,36);
    three_angle_range=nchoosek(angle_range,2);
    normdist_range=linspace(0.1,1,10);
    append_range=normdist_range(1).*ones(length(three_angle_range),1);
    final_range=horzcat(append_range,three_angle_range);
    for i=2:length(normdist_range)
        append_range=normdist_range(i).*ones(length(three_angle_range),1);
        final_range=vertcat(final_range,horzcat(append_range,three_angle_range));
    end
    f1 = figure('visible','on');
    final_range2=ml_mappos(final_range);
    y_a=1./(1+exp(-final_range2(:,2:6)*spatial.beta));
    [x_1,y_1,z_1]=sph2cart(final_range(:,2),final_range(:,3),final_range(:,1));
    scatter3(x_1,y_1,z_1,'filled','CData',y_a);
    title('Spatial Distribution of Model  fitted to a Logistic Regression Model');
    colorbar;
    saveas( f1, 'spatial_model_fit.png', 'png');
    I = imread( 'spatial_model_fit.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'spatial_model_fit_thumbnail.png' );
    img2html(fileID,'spatial_model_fit.png','spatial_model_fit.png','Spatial Distribution fitted to a Logistic Regression Model');


end