function [Overlay] = generate_overlay_from_gfp_dic(GFP, DIC)
% Generate overlay image from GFP and DIC image


img_num = numel(GFP);
Overlay = cell(img_num, 1);

for i = 1 : img_num
    cur_overlay = repmat(double(DIC{i}), [1, 1, 3]) .* reshape([1, 0.7, 1], [1, 1, 3]);
    cur_overlay(:, :, 2) = cur_overlay(:, :, 2) + double(GFP{i}) * 0.4;
    
    Overlay{i} = uint8(cur_overlay);
end

end
