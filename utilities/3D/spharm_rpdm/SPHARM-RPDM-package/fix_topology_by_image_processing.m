function [bim_out, vertices, faces] = fix_topology_by_image_processing(bim)
% fix topology to make sure Euler character is 2 for the mesh
% use isosurface for judgement, which is faster.

% 1/11/2021 R.F. Murphy trap cases where eroded image goes to zero
% 4/24/2023 R.F. Murphy trap more such cases and fix message logic

origin = [0, 0, 0];
vxsize = [1, 1, 1];
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
disp(['fix_topology_by_image_processing 1st try; voxels,vertices,faces:',...
    num2str(sum(bim(:))),',',num2str(size(vertices,1)),',',...
    num2str(size(faces,1))])
bim_out = bim;

addpath(genpath('image_geometric_measures'));
if size(vertices, 1) - size(faces, 1) ~= 2
    se_size = 0; bim_1 = bim;
    [chi, labels] = imEuler3d(bim_1);
    while chi ~= 1 && se_size <= 20
        se_size = se_size + 1;
        se = strel('disk', se_size);
        bim_1 = imdilate(bim, se);
        bim_1 = imerode(bim_1, se);
        % check whether everything was eroded meaning approach didn't work
        if sum(bim_1(:))<=0
            bim_1 = bim;
            break
        end
        [chi, labels] = imEuler3d(bim_1);
    end
    % only use the connected component that fixes the holes. 
    if chi == 1
        bim_diff = bim_1 - bim;
        CC = bwconncomp(bim_diff > 0, 6);
        bim_2 = bim_1;
        for i = 1 : CC.NumObjects
            cur_pixel_inds = CC.PixelIdxList{i};
            bim_2(cur_pixel_inds) = false;
            [chi, labels] = imEuler3d(bim_2);
            if chi ~= 1
                bim_2(cur_pixel_inds) = true;
            end     
        end
        [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
        bim_out = bim_2;
    end
end
% if still cannot fix the topology, then use gen_surf_data as criterion
% and use imopen to remove some small voxels
if size(vertices, 1) - size(faces, 1) ~= 2
    disp(['fix_topology_by_image_processing 2nd try; voxels,vertices,faces:',...
        num2str(sum(bim_out(:))),',',num2str(size(vertices,1)),',',...
        num2str(size(faces,1))])
    se_size = 0; bim_1 = bim;
    [vertices, faces] =  gen_surf_data(bim_1,origin,vxsize);
    while  size(vertices, 1) - size(faces, 1) ~= 2 && se_size <= 20
        se_size = se_size + 1;
        se = strel('disk', se_size);
        bim_1_prev = bim_1;
        bim_1 = imerode(bim, se);
        % check whether everything was eroded meaning approach didn't work,
        % in which case keep the previous iteration
        if sum(bim_1(:))<=0
            bim_1 = bim_1_prev;
            break
        end
        bim_1 = imdilate(bim_1, se);
        CC = bwconncomp(bim_1 > 0);
        [~, max_ind] = max(cellfun(@numel, CC.PixelIdxList));
        bim_1 = false(size(bim));
        bim_1(CC.PixelIdxList{max_ind}) = true;
        [vertices, faces] =  gen_surf_data(bim_1,origin,vxsize);
        % this shouldn't happen
        if isempty(vertices)
            bim_1 = bim;
            break
        end
    end
    bim_diff = bim - bim_1;
    CC = bwconncomp(bim_diff > 0, 6);
    bim_2 = bim_1;
    for i = 1 : CC.NumObjects
        cur_pixel_inds = CC.PixelIdxList{i};
        bim_2(cur_pixel_inds) = true;
        [chi, labels] = imEuler3d(bim_2);
        if chi ~= 1
            bim_2(cur_pixel_inds) = false;
        end     
    end
    
    [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
    bim_out = bim_2;  
    % if still not working, then use the Euler number as criterion, as
    % generation of mesh is slow, we use a group update, that is, we first
    % put N CCs together, if not work, then evaluate them and remove bad ones.
    if size(vertices, 1) - size(faces, 1) ~= 2
        disp(['fix_topology_by_image_processing 3rd try; voxels,vertices,faces:',...
            num2str(sum(bim_out(:))),',',num2str(size(vertices,1)),',',...
            num2str(size(faces,1))])
        bim_2 = bim_1;
        per_partition = 20;
        if CC.NumObjects < 10
            per_partition = 1;
        end
        N_partition = ceil(CC.NumObjects / per_partition);            
%         pixel_nums = cellfun(@numel, CC.PixelIdxList);
%         [~, pixel_inds] = sort(pixel_nums, 'descend');
%         CC.PixelIdxList = CC.PixelIdxList(pixel_inds);
        for i = 1 : N_partition
            ind_i = (i-1) * per_partition + 1 : min(CC.NumObjects, i * per_partition);
            cur_pixel_cells = CC.PixelIdxList(ind_i);
            cur_pixel_inds = cat(1, cur_pixel_cells{:});
            bim_2(cur_pixel_inds) = true;
            [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
            if size(vertices, 1) - size(faces, 1) ~= 2
                for j = 1 : numel(ind_i)
                    cur_pixel_inds_j = CC.PixelIdxList{ind_i(j)};
                    bim_2(cur_pixel_inds_j) = ~true;
                    [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
                    if size(vertices, 1) - size(faces, 1) == 2
                        break;
                    else
                        bim_2(cur_pixel_inds_j) = true;
                    end  
                    % if at the end, the topology still not valid, remove
                    % all voxel in this partition.
                    if j == numel(ind_i)
                        bim_2(cur_pixel_inds) = ~true;
                    end                        
                end
            end     
        end
        [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
        bim_out = bim_2;
    end
end

CC = bwconncomp(bim_out > 0, 6);
% if still not fixed, then first close then open to see if it works
if size(vertices, 1) - size(faces, 1) ~= 2 || CC.NumObjects ~= 1
    disp(['fix_topology_by_image_processing 4th try; voxels,vertices,faces:',...
        num2str(sum(bim_out(:))),',',num2str(size(vertices,1)),',',...
        num2str(size(faces,1))])
    se_size = 0; bim_1 = bim;
    [vertices, faces] =  gen_surf_data(bim_1,origin,vxsize);
    while  size(vertices, 1) - size(faces, 1) ~= 2 && se_size <= 20
        se_size = se_size + 1;
        se = strel('disk', se_size);
        bim_1_prev = bim_1;
        bim_1 = imdilate(bim, se);
        bim_1 = imerode(bim_1, se);
        bim_1 = imerode(bim_1, se);
        bim_1 = imdilate(bim_1, se);
        if sum(bim_1(:))<=0
            bim_1 = bim_1_prev;
            break
        end
        [vertices, faces] =  gen_surf_data(bim_1,origin,vxsize);
        % this shouldn't happen
        if isempty(vertices)
            bim_1 = bim;
            break
        end

    end
    bim_diff = bim - bim_1;
    CC = bwconncomp(bim_diff > 0, 6);
    bim_2 = bim_1;
    for i = 1 : CC.NumObjects
        cur_pixel_inds = CC.PixelIdxList{i};
        bim_2(cur_pixel_inds) = true;
        [chi, labels] = imEuler3d(bim_2);
        if chi ~= 1
            bim_2(cur_pixel_inds) = false;
        end     
    end
    [vertices, faces] =  gen_surf_data(bim_2,origin,vxsize);
    bim_out = bim_2;    
    if size(vertices, 1) - size(faces, 1) ~= 2
        bim_out = bim_1;
    end    
end
% check if everything failed and return the input image
if size(vertices, 1) - size(faces, 1) ~= 2
    disp('fix_topology_by_image_processing failed; returning original image')
    bim_out = bim;
end
[vertices, faces] =  gen_surf_data(bim_out,origin,vxsize);

end