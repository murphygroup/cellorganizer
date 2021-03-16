function [ip] = image_inner_product(ax, ay, az, bx, by, bz)
  ip = ax .* bx + ay .* by + az .* bz;
end
