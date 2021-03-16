function [I_out]= del2_3d(I_input)
% the discrete laplacian operation for a 3d image. This function is faster
% than the built-in one. 
% Authore: xruan 
% Date

    if ndims(I_input) ~= 3
        error('the input image must be in 3d');
    end

    [l, m, n] = size(I_input);

    f = I_input;

    gr = zeros(l, m, n);

    gr(2 : l - 1, :, :) = (diff(f(2 : l, :, :), 1, 1) - diff(f(1 : l - 1, :, :), 1, 1)) ./ 2;

    gr(1, :, :) = gr(2, :, :) .* 2 - gr(3, :, :);
    gr(l, :, :) = gr(l - 1, :, :) * 2 - gr(l - 2, :, :);

    gc = zeros(l, m, n);

    gc(:, 2 : m - 1, :) = (diff(f(:, 2 : m, :), 1, 2) - diff(f(:, 1 : m - 1, :), 1, 2)) ./ 2;

    gc(:, 1, :) = gc(:, 2, :) .* 2 - gc(:, 3, :);
    gc(:, m, :) = gc(:, m - 1, :) * 2 - gc(:, m - 2, :);


    gs = zeros(l, m, n);

    gs(:, :, 2 : n - 1) = (diff(f(:, :, 2 : n), 1, 3) - diff(f(:, :, 1 : n - 1), 1, 3)) ./ 2;

    gs(:, :, 1) = gs(:, :, 2) .* 2 - gs(:, :, 3);
    gs(:, :, n) = gs(:, :, n - 1) * 2 - gs(:, :, n - 2);

    I_out = (gr + gc + gs) / 3;


end