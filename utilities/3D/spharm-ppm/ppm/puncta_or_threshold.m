function [processed_imgs,puncta]=puncta_or_threshold(imgs,puncta_or_threshold_flag)
    
% %%% run such command in matlab
% imgs = cell(1,1);
% imgs{1} = bfOpen3DVolume('/Users/piepi/CellOrganizer/Josh/test2.tif');
% channels = cell(5,length(imgs));
% i = 1;
% z = size(imgs{i}{1}{1},3)/5;
% for j = 1:5
%     % j*z-z+1
%     % j*z
%     channels{j,i} = imgs{i}{1}{1}(:,:,j*z-z+1:j*z); %%% 这里可以改一下成preprocess_2的 ？？
% end
% imgs = channels;
% prot_imgs = imgs;
% seg_channels=prot_imgs;
% chan_regs=seg_channels;

%% below are script belong to this function
z_scale = 0.4/0.161427354511474;
improt=imgs; %% a 3D channel that has multiple z planes

%% threshold image 
f1=1.5;
f2=5;
% function [fg,bg, fg_denoise] = segmentProt(improt,f1, f2)
fg_denoise = zeros(size(improt));
for z = 1:size(improt,3)
    curslice = double(improt(:,:,z));

    f = fspecial('gaussian', 15, f1);
    improt_f = imfilter(curslice,f);
    f = fspecial('gaussian', 15, f2);
    im_f = imfilter(curslice, f);

    improt_fg = improt_f - im_f;
    improt_fg(improt_fg < 0) = 0;

    improt_fg_denoise = improt_fg;
    improt_fg_denoise(improt_fg_denoise < ml_rcthreshold(uint8(improt_fg_denoise))) = 0;
    improt_fg_denoise = improt_fg_denoise.* bwmorph(improt_fg_denoise, 'clean', inf);

    fg_denoise(:,:,z) = improt_fg_denoise;
end % End of function
pdn=fg_denoise;

%find obj, threshold and remain big ones, and remove small ones in each plane
disp(['ml_findobjs...'])
[p1ob] = ml_findobjs(pdn);
p1obfinal=[];
k=1;
for i=1:length(p1ob)
    if (size(p1ob{i},1)>=4)
        p1obfinal{k}=p1ob{i};
        k=k+1;
    else
        pdn(p1ob{i}(:,1),p1ob{i}(:,2),p1ob{i}(:,3)) = 0;
    end
end
%pdn here is obj in each plane after threshhold
if puncta_or_threshold_flag==0
    processed_imgs=pdn;
    puncta=[];
    disp(['image threshhold complete'])
elseif puncta_or_threshold_flag==1
    %find centers
    p1labels=bwlabeln(pdn,26);
    objnum=max(p1labels(:));
    p1centers=[];
    for k = 1:objnum
        k
        objnum
        p1posxyz=find(p1labels==k);
        [p1X,p1Y,p1Z] = ind2sub(size(pdn),p1posxyz);
        v=pdn(p1labels == k);
        pos = find(v==max(v));
        % p1centers = [p1centers;[p1X(pos) p1Y(pos) p1Z(pos)*z_scale]];
        p1centers = [p1centers;[p1X(pos) p1Y(pos) p1Z(pos)]];        
    end
    puncta=p1centers;
    % get puncta image
    processed_imgs=zeros(size(imgs));
    for i=1:size(puncta,1)
        processed_imgs(puncta(i,1),puncta(i,2),puncta(i,3))=1;
    end
    disp(['puncta finding complete'])
else
    warning('Wrong puncta_or_threshold_flag. puncta_or_threshold_flag should be 0 or 1 as int.');
    return;
end

end