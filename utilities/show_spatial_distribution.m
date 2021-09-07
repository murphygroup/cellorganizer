function [ output_args ] = show_spatial_distribution( model,fileID)
%Author : Serena Abraham



    f = figure('visible','on');
    spatial=model;
    [x1,y1,z1]=sph2cart(spatial.anglestheta,spatial.anglesphi,spatial.normdists);
    [X,Y,Z] = sphere;
    [az,phi,r]=cart2sph(X(:),Y(:),Z(:));
    new_az=[];new_phi=[];new_r=[];
    for i=1:length(az)
        if and(phi(i)<=max(spatial.anglesphi),az(i)<=max(spatial.anglestheta))
            if and(phi(i)>=min(spatial.anglesphi),az(i)>=min(spatial.anglestheta))
                new_az(end+1)=az(i);
                new_phi(end+1)=phi(i);
                new_r(end+1)=r(i);
            end
        end
    end
    [new_x,new_y,new_z]=sph2cart(new_az,new_phi,new_r);
    
    scatter3(new_x,new_y,new_z,'filled','MarkerFaceColor',[0.7 0.7 0.7]);
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
    %angle_range=linspace(0,6.28319,36);
    %three_angle_range=nchoosek(angle_range,2);
    %normdist_range=linspace(0.1,1,10);
    %append_range=normdist_range(1).*ones(length(three_angle_range),1);
    %final_range=horzcat(append_range,three_angle_range);
    %for i=2:length(normdist_range)
    %    append_range=normdist_range(i).*ones(length(three_angle_range),1);
    %    final_range=vertcat(final_range,horzcat(append_range,three_angle_range));
    %end
    %f1 = figure('visible','off');
    [r1,r2,r3]=cart2sph(new_x(:),new_y(:),new_z(:));
    %final_range2=ml_mappos(final_range);
    final_range2=ml_mappos(horzcat(r3,r2,r1));
    %y_a=1./(1+exp(-final_range2(:,2:6)*spatial.beta));
    %[x_1,y_1,z_1]=sph2cart(final_range(:,2),final_range(:,3),final_range(:,1));
    %scatter3(x_1,y_1,z_1,'filled','CData',y_a);
    %title('Relative Probability of a Spatial Positions using a Logistic Regression Model');
    %colorbar;
    %saveas( f1, 'relative_probability.png', 'png');
    %I = imread( 'relative_probability.png' );
    %I = imresize( I, 0.50 );
    %imwrite( I, 'relative_probability_thumbnail.png' );
    %img2html(fileID,'relative_probability.png','relative_probability_thumbnail.png','Relative Probability of a Spatial Positions');
    
    f2 = figure('visible','on');
    y_b = mvksdensity(spatial.mappos_x(:,2:6),final_range2(:,2:6),'Bandwidth',0.28,'Kernel','normal');
    scatter3(new_x,new_y,new_z,'filled','CData',y_b);
    title('Spatial Probability Distribution of Model  fitted to a Gaussian Kernel');
    colorbar;
    saveas( f2, 'spatial_model_fit.png', 'png');
    I = imread( 'spatial_model_fit.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'spatial_model_fit_thumbnail.png' );
    img2html(fileID,'spatial_model_fit.png','spatial_model_fit_thumbnail.png','Spatial Probability Distribution fitted to a Gaussian Kernel');



end