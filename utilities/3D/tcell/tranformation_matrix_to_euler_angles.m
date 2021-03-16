function [theta_1, theta_2, theta_3, scale_factor] = tranformation_matrix_to_euler_angles(M, order)
% the function extract euler angles and scale factor from transformation
% matrix. theta_1, 2, 3 represents angle of x, y, z -axis, with equation
% R_x(theta_1) R_y(theta_2) R_z(theta_3)
% https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf
% 
% Author: Xiongtao Ruan
% Date: Aug. 08, 2016
%
% Sept.1, update another order for Euler angles, 'ZXY'

if nargin < 2
    order = 'ZXY';
end

M_3d = M(1 : 3, 1 : 3);
scale_factor = det(M_3d) ^ (1 / 3);

M_r = M_3d / scale_factor;

eps_thrsh = 5 * eps;

switch order 
    case 'ZYX'
        theta_1 = atan2(M_r(2, 3), M_r(3, 3));
        c_2 = sqrt(M_r(1, 1) .^ 2 + M_r(1, 2) .^ 2);
        theta_2 = atan2(-M_r(1, 3), c_2);
        s_1 = sin(theta_1);
        c_1 = cos(theta_1);
        theta_3 = atan2(s_1 * M_r(3, 1) - c_1 * M_r(2, 1), c_1 * M_r(2, 2) - s_1 * M_r(3, 2));
    case 'ZXY'
        
        % threshold the floating error
        M_r(abs(M_r) < eps_thrsh) = 0;
        
        if M_r(2, 1) == 0 && M_r(2, 2) == 0
            if 1 - abs(M_r(2, 3)) > eps_thrsh
                error('invalid rotation matrix!');
            end
            theta_1 = atan2(M_r(2, 3), 0);
            theta_2 = 0;
            theta_3 = atan2(M_r(1, 2), M_r(1, 1));
        elseif M_r(2, 1) == 0
            % first try to only use theta_3 = 0, and see if it is correct
            theta_3 = 0;
            theta_2 = atan2(M_r(3, 1), M_r(1, 1));
            theta_1 = atan2(M_r(2, 3), M_r(2, 2));
        else
            theta_3 = atan2(-M_r(2, 1), M_r(2, 2));
            c_1 = -M_r(2, 1) ./ sin(theta_3);
            theta_1 = atan2(M_r(2, 3), c_1);
            s_2 = - M_r(1, 3) / c_1;
            c_2 = M_r(3, 3) / c_1;
            theta_2 = atan2(s_2, c_2);            
        end
        
        if false 
            theta_3 = atan2(-M_r(2, 1), M_r(2, 2));

            c_1 = sqrt(M_r(1, 3) .^ 2 + M_r(3, 3) .^ 2);
            if abs(sin(theta_3)) > eps_thrsh 
                c_1 = sign(-M_r(2, 1)) * sign(sin(theta_3)) * c_1;
            else
                c_1 = sign(M_r(2, 2)) * sign(cos(theta_3)) * c_1;
            end
            theta_1 = atan2(M_r(2, 3), c_1);
            if theta_1 == 0 && theta_3 == 0
                theta_2 = atan2(-M_r(1, 3), M_r(3, 3));
            else
                s_3 = sin(theta_3);
                c_3 = cos(theta_3);
                theta_2 = atan2(c_3 * M_r(1, 2) - s_3 * M_r(1, 1), s_3 * M_r(3, 1) - c_3 * M_r(3, 2));    
            end
        end
        
end
           

end


