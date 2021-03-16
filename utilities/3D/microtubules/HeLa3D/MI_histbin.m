function [Histogram] = MI_histbin(image)
% Histogram features .. 32 bins

image = uint8(image);
image = image(image~=0);

Histogram = histc(image(:),1:255)';
Histogram = [Histogram,hist(image(:),2:3:255),hist(image(:),4:7:255)];

end % End of function
