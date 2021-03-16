function features = MI_feat_extract_truncate(input_image,featArray)

% Truncate input image to 255
input_image = uint8(input_image);

% STEP 4 - Extract texture features


direction = int32([0 0 1;
0 1 0; 0 1 1; 0 1 -1;
1 0 0; 1 0 1; 1 0 -1;
1 1 0; 1 1 1; 1 1 -1;
1 -1 0; 1 -1 1; 1 -1 -1]);

[f,n] = ml_3Dcoocmat(input_image,direction); f(1,1) = 0;
nb = int32(length(n));
f = f./sum(f(:));
[features,featuresnames] = ml_haralicktexture(single(f),n,nb);

features = features(featArray)';

end
