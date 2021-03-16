function save_spharm_obj(obj_img,min_obj_size,cellName)

objs = bwconncomp(obj_img);

% BW(objs.PixelIdxList{idx(i)}) = 0;
save_dir = [pwd filesep 'spharm_input' ];
[status, msg, msgID] = mkdir(save_dir);
for idx = 1:length(objs.PixelIdxList)
    tmp_obj = zeros(size(obj_img));
    if length(objs.PixelIdxList{idx})<min_obj_size
        continue;
    end
    tmp_obj(objs.PixelIdxList{idx})=1;
    [x,y,z] = ind2sub(size(obj_img),objs.PixelIdxList{idx}); % find the small box to extract image
    sz = size(tmp_obj);
    [x_down,x_up] = get_bound(x,sz(1));
    [y_down,y_up] = get_bound(y,sz(2));
    [z_down,z_up] = get_bound(z,sz(3));
    sf_img = tmp_obj(x_down:x_up,y_down:y_up,z_down:z_up);
%     disp(size(sf_img));
%     disp(numel(objs.PixelIdxList{idx}));
    sfp = [save_dir filesep cellName '_obj_' int2str(idx) '.tif'];
    save_3D_tif(sf_img,sfp);
    
end
end

% imwrite(im1,'myMultipageFile.tif')
% imwrite(im2,'myMultipageFile.tif','WriteMode','append')

% help function to save 3D tif file from a matrix
function save_3D_tif(img,sfp)
sz = size(img);

for i=1:sz(3)
  tiff = img(:, :, i);
  if i==1
      imwrite(tiff,sfp);
  end
  imwrite(tiff,sfp,'WriteMode', 'append')
end
end

function [x_down,x_up]=get_bound(x,max_x)
x_down = max(1,min(x)-1);
x_up = min(max(x)+1,max_x);
end