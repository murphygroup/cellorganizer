function [ PSF ] = get_mt_psf()

rad = 240;

bb = imread(['SVI_PSF_XZ_' num2str(rad) '.jpg']);

bb(385:392,1:16) = 0;
b = imresize(bb,[size(bb,1)/9,size(bb,2)/3],'nearest');

b = imresize(b,[size(b,1)*(0.148/0.2), size(b,2)]);

% [cx,cy,sx,sy,PeakOD] = gaussian2d(double(a),1);
midpt = round(size(b,2)/2);
K = 1;

for idx = 1:size(b,1)

        m = double(b(idx,:));
        if sum(m(:))>0
                m = m./sum(m(:));
                vv1 = repmat(m,size(m,2),1);
                vv2 = repmat(m',1,size(m,2));

                vv3 = vv1.*vv2;
                vv3 = vv3./sum(vv3(:));
                vv3 = (double(max(b(idx,:)))/max(vv3(:)))*vv3;

                result_all(:,:,K) = vv3;
                K = K + 1;
        end
end

result_all(1,:,:) = [];
result_all(:,1,:) = [];

result_all = ml_3dXYresize(result_all,(0.049/0.2));

result_all = uint8(result_all);

[I,J,K] = ind2sub(size(result_all),find(result_all(:) == max(result_all(:))));

result_all(:,:,[1:(K-6),(K+6):end]) = [];

PSF = double(result_all);
PSF = PSF./sum(PSF(:));

end

