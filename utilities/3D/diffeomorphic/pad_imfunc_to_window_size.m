function img = pad_imfunc_to_window_size(imgnum, imfunc, window_size )
%takes a function that returns a matrix of arbitrary size, and pads the x
%and y dimensions to be divisible by the window_size

img = imfunc(imgnum);

if ~isempty(img)

    imsize = size(img);
    imsize = imsize(1:2);

    targetsize = ceil(max(imsize)/window_size) * window_size;

    padding = (targetsize - imsize)/2;

    img = padarray(img, [floor(padding) 0]);
    img = padarray(img, [abs(floor(padding) - ceil(padding)) 0], 'pre');

end

end
