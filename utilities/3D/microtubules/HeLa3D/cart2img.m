function image_coord = cart2img(newpoint, imgcent_coordinate, XYZres)

image_coord = [(-round((1/XYZres)*newpoint(2)) + imgcent_coordinate(1)); ...
                (round((1/XYZres)*newpoint(1)) + imgcent_coordinate(2));...
                (-round((1/XYZres)*newpoint(3)) + imgcent_coordinate(3))];
            
end
