function [spharm_descriptors] = convert_3d_preshape_to_spharm(preshapes, scales, centers)
% convert preshape to spharm descriptor with original centers and scales

ln = size(preshapes, 1);
preshapes_1 = preshapes .* permute(scales, [2, 3, 1]);
preshapes_2 = reshape(preshapes_1, ln / 2, 2, size(preshapes_1, 2), size(preshapes_1, 3));
preshapes_3 = squeeze(preshapes_2(:, 1, :, :) + 1i .* preshapes_2(:, 2, :, :));
preshapes_4 = preshapes_3;
preshapes_4 = cat(1, centers, preshapes_4);
spharm_descriptors = preshapes_4;

end

