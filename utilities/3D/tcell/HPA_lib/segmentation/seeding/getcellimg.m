function img = getcellimg( img)

img = imfilter(img,[3 3]);

[c b] = imhist(img);
[a ind] = max(c);

st = ceil(std(double(img(img>min(img(:))))));

% cell = zeros(size(img));
% for V = 0:st
%     cell = cell+(img>b(ind)+V);
% end
% cell = uint8(cell);

img(img>st) = st+1;
