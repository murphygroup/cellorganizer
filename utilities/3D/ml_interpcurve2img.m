function sliceimg = ml_interpcurve2img( imgsize, x, y )

    sliceimg = zeros(imgsize);
    any_points_outside_image = false;
    % interpolate using 100 points for each line segment
    for t = 1:length(x)-1
        rpts = round(linspace(y(t),y(t+1),100));
        cpts = round(linspace(x(t),x(t+1),100));

        % T. Buck 2019-01-29 Prevent "Out of range subscript" error in sub2ind without issuing a warning for every t or skipping any pixels that should be in the image
        %{
        points_inside_image = (rpts >= 0.5) & (cpts >= 0.5) & (rpts < imgsize(1)) & (cpts < imgsize(2));
        any_points_outside_image = any_points_outside_image || ~all(points_inside_image(:));
        rpts = rpts(points_inside_image);
        cpts = cpts(points_inside_image);
        %}
        % Prevent "Out of range subscript." error but keep filling functionality
        rpts = max(rpts, 1);
        rpts = min(rpts, imgsize(1));
        cpts = max(cpts, 1);
        cpts = min(cpts, imgsize(2));

        index = sub2ind(imgsize,rpts,cpts);
        sliceimg(index) = 255;
    end
    if any_points_outside_image
        warning('CellOrganizer:ml_interpcurve2img:outOfRange', 'CellOrganizer: Out of range subscripts when adding pixels to image');
    end

end

