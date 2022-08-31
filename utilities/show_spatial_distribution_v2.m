function KL_div = show_spatial_distribution_v2( model,fileID)
%Author : Serena Abraham, Paul Sun


    KL_div = 0;
    f = figure('visible','on');
    spatial=model;
    if length(model) > 1
        normdists=horzcat(model{1}.normdists, model{2}.normdists);
        anglestheta=horzcat(model{1}.anglestheta, model{2}.anglestheta);
        anglesphi=horzcat(model{1}.anglesphi, model{2}.anglesphi);
        [x1,y1,z1]=sph2cart(anglestheta, anglesphi, normdists);
    else
    [x1,y1,z1]=sph2cart(spatial.anglestheta,spatial.anglesphi,spatial.normdists);
    end
    [X,Y,Z] = sphere;
    [az,phi,r]=cart2sph(X(:),Y(:),Z(:));
%     new_az=[];new_phi=[];new_r=[];
%     for i=1:length(az)
%         if and(phi(i)<=max(spatial.anglesphi),az(i)<=max(spatial.anglestheta))
%             if and(phi(i)>=min(spatial.anglesphi),az(i)>=min(spatial.anglestheta))
%                 new_az(end+1)=az(i);
%                 new_phi(end+1)=phi(i);
%                 new_r(end+1)=r(i);
%             end
%         end
%     end
    new_az=az;
    new_phi=phi;
	new_r=r;
    [new_x,new_y,new_z]=sph2cart(new_az,new_phi,new_r);
    
    if length(model) > 1
        d1 = length(model{1}.normdists);
        d2 = length(model{2}.normdists);
        scatter3(x1(1:d1),y1(1:d1),z1(1:d1),'filled','MarkerFaceColor', 'blue');
        hold on
        scatter3(x1(d1+1:end),y1(d1+1:end),z1(d1+1:end),'filled','MarkerFaceColor', 'cyan');
        hold on
        scatter3(new_x,new_y,new_z,'filled','MarkerFaceColor', [0.7, 0.7, 0.7]);
        legend('model 1','model 2','cell boundary');
    else
        scatter3(x1,y1,z1,'filled','b');
        hold on
        scatter3(new_x,new_y,new_z,'filled','MarkerFaceColor', [0.7, 0.7, 0.7]);
        legend('model', 'cell boundary');
    end
    hold on
    
   
    
    title('Spatial Distribution of objects');
    hold off;
    saveas( f, 'show_spatial_distribution.png', 'png');
    I = imread( 'show_spatial_distribution.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_spatial_distribution_thumbnail.png' );
    img2html(fileID,'show_spatial_distribution.png','show_spatial_distribution_thumbnail.png','Spatial Distribution across both models');
    
    az_ = [];
    phi_ = [];
    r_ = [];
    for i=1:10
        az_ = [az_ az'];
        phi_ = [phi_ phi'];
        r_ = [r_ 0.1*i*r'];
    end
    
%     [r1,r2,r3]=cart2sph(az_',phi_',r_');
    [new_x,new_y,new_z]=sph2cart(az_',phi_',r_');
    %final_range2=ml_mappos(final_range);
    final_range2=ml_mappos(horzcat(new_x,new_y,new_z));
   
    f2 = figure('visible','on');
%     y_b = mvksdensity(spatial.mappos_x(:,2:6),final_range2(:,2:6),'Bandwidth',0.28,'Kernel','normal');
    if length(model) > 1
        y_b1 = mvksdensity([x1(1:d1)' y1(1:d1)' z1(1:d1)'],[new_x new_y new_z],'Bandwidth',0.28,'Kernel','normal');
        y_b2 = mvksdensity([x1(d1+1:end)' y1(d1+1:end)' z1(d1+1:end)'],[new_x new_y new_z],'Bandwidth',0.28,'Kernel','normal');
        scatter3(new_x,new_y,new_z, 20, 'filled','CData',y_b2);
        title('Spatial Probability Distribution of Model 2 fitted to a Gaussian Kernel');
        colorbar;
        f3 = figure('visible','on');
        scatter3(new_x,new_y,new_z, 20, 'filled','CData',abs(y_b2-y_b1));
        title('Difference in Spatial Probability Distribution between 2 Models fitted to a Gaussian Kernel');
        colorbar;
        
        
    else
        y_b = mvksdensity([x1' y1' z1'],[new_x new_y new_z],'Bandwidth',0.28,'Kernel','normal');
        scatter3(new_x,new_y,new_z, 20, 'filled','CData',y_b);
        title('Spatial Probability Distribution of Model fitted to a Gaussian Kernel');
        colorbar;
    end
    saveas( f2, 'spatial_model_fit.png', 'png');
    I = imread( 'spatial_model_fit.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'spatial_model_fit_thumbnail.png' );
    img2html(fileID,'spatial_model_fit.png','spatial_model_fit_thumbnail.png','Spatial Probability Distribution fitted to a Gaussian Kernel');

    if length(model) > 1
        saveas( f3, 'spatial_model_diff.png', 'png');
        I = imread( 'spatial_model_diff.png' );
        I = imresize( I, 0.50 );
        imwrite( I, 'spatial_model_diff_thumbnail.png' );
        img2html(fileID,'spatial_model_diff.png','spatial_model_diff_thumbnail.png','Difference in Spatial Probability Distribution between 2 Models fitted to a Gaussian Kernel');
        y_b1 = y_b1 + 1e-5;
        y_b2 = y_b2 + 1e-5;
        y_b1 = y_b1/sum(y_b1);
        y_b2 = y_b2/sum(y_b2);
        KL_div = sum(y_b1.*log(y_b1./y_b2));
        text2html(fileID, sprintf('KL divergence between spatial distributions of model1 and model2: %.4f;\n', KL_div));
        
    end
end