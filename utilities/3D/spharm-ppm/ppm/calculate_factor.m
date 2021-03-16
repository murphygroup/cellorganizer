function factor_matrix=calculate_factor(img_path,model_channel_img,other_channel_imgs,puncta_or_threshold_flags,aspect)

    [~,mc_puncta]=puncta_or_threshold(model_channel_img,1);
    save([img_path '__mc_puncta.mat'],'mc_puncta');
    % mc_puncta=load([img_path '__mc_puncta.mat']);
    % mc_puncta=mc_puncta.mc_puncta;

    %calculate distance transformation matrix for other_channel_imgs
    dist_trans_matrices={};
    for j=1:size(other_channel_imgs,2)
        [oc_processed_imgs,~]=puncta_or_threshold(other_channel_imgs{j},puncta_or_threshold_flags(j));
        save([img_path '__oc_processed_imgs' int2str(j) '.mat'],'oc_processed_imgs');
        dist_trans_matrices{j}=bwdistsc(oc_processed_imgs,aspect); %%% remember to cite this in paper
        temp_dist_trans_matrices=dist_trans_matrices{j};
        save([img_path '__dist_trans_matrices' int2str(j) '.mat'],'temp_dist_trans_matrices');
    end

    %calculate factor matrix
    factor_matrix=zeros(size(mc_puncta,1),size(other_channel_imgs,2));
    for i=1:size(mc_puncta,1)
        for j=1:size(dist_trans_matrices,2)
            factor_matrix(i,j)=dist_trans_matrices{j}(mc_puncta(i,1),mc_puncta(i,2),mc_puncta(i,3));
        end
    end
    save([img_path '__factor_matrix.mat'],'factor_matrix');
end