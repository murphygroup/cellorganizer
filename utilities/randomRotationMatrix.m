function rotation_matrix = randomRotationMatrix(dimensionality)
% RANDOMROTATIONMATRIX Generates a random rotation matrix.

% Taraz Buck 2019-04-09

switch dimensionality
    case 3
        % Arvo 1992, "Fast Random Rotation Matrices"
        x1 = rand * pi;
        R = [cos(x1) sin(x1) 0; -sin(x1) cos(x1) 0; 0 0 1];
        x2 = rand;
        x3 = rand;
        v = [cos(2 * pi * x2) * sqrt(x3); sin(2 * pi * x2) * sqrt(x3); sqrt(1-x3)];
        H = eye(3) - 2 * v * v';
        M = -H * R;
        rotation_matrix = M;
    otherwise
        error('Dimensionality %i not implemented');
end

