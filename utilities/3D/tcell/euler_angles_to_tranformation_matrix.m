function [M] = euler_angles_to_tranformation_matrix(theta_1, theta_2, theta_3, scale_factor, order)
% the function use eular angles and scale factor to compute transformation
% matrix. theta_1, 2, 3 represents angle of x, y, z -axis, with equation
% R_x(theta_1) R_y(theta_2) R_z(theta_3) R_scale for ZYX
% https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf
% 
% Author: Xiongtao Ruan
% Date: Aug. 08, 2016

if nargin < 4
    scale_factor = 1;
end
if nargin < 5
    order = 'ZXY'
end

C = cos([theta_1, theta_2, theta_3]);
S = sin([theta_1, theta_2, theta_3]);
switch order
    case 'ZYX'
        M = [C(2) * C(3),  C(2) * S(3),  - S(2);
             S(1) * S(2) * C(3) - C(1) * S(3),  S(1) * S(2) * S(3) + C(1) * C(3), S(1) * C(2);
             C(1) * S(2) * C(3) + S(1) * S(3), C(1) * S(2) * S(3) - S(1) * C(3), C(1) * C(2)];
    case 'ZXY'
        M = [C(2) * C(3) - S(1) * S(2) * S(3),  C(2) * S(3) + C(3) * S(1) * S(2),  -C(1) * S(2);
             -C(1) * S(3), C(1) * C(3), S(1);
             C(3) * S(2) + C(2) * S(1) * S(3), S(2) * S(3) - C(2) * C(3) * S(1), C(1) * C(2)];
end
        
M = M * scale_factor;

end
