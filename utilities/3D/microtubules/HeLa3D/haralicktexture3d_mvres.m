function features = haralicktexture3d_mvres(input_image,featArray)

% Image res in X and Y is the same resolution. Only Z is different resolution.

% This code divides the 13 directions into 2.

% Truncate input image to 255

input_image = uint8(input_image);

% STEP 4 - Extract texture features


direction1 = int32([0 0 1;
0 1 0; 0 1 1; 0 1 -1]);

direction2 = int32([1 0 0;
1 0 1; 1 0 -1;
1 1 0; 1 1 1; 1 1 -1;
1 -1 0; 1 -1 1; 1 -1 -1]);

[f1,n1] = ml_3Dcoocmat(input_image,direction1); f1(1,1) = 0;
nb1 = int32(length(n1));
f1 = f1./sum(f1(:));
[features1,featuresnames1] = ml_haralicktexture(single(f1),n1,nb1);

[f2,n2] = ml_3Dcoocmat(input_image,direction2); f2(1,1) = 0;
nb2 = int32(length(n2));
f2 = f2./sum(f2(:));
[features2,featuresnames2] = ml_haralicktexture(single(f2),n2,nb2);

features = [features1(featArray)',features2(featArray)'];

end
