% test symbolic computing of gradient and hessian. 

x = zeros(3, 4);

x = sym('x%d%d', [3, 4], 'real');

a = x(:, 1);
b = x(:, 2);
c = x(:, 3);
d = x(:, 4);

ab = dot(a, b);
ac = dot(a, c);
ad = dot(a, d);
bc = dot(b, c);
bd = dot(b, d);
cd = dot(c, d);

Ca = bd - ad .* ab;
Cb = ac - ab .* bc;
Cc = bd - bc .* cd;
Cd = ac - cd .* ad;

spat(1) = det([d, a, b]);
spat(2) = det([a, b, c]);
spat(3) = det([b, c, d]);
spat(4) = det([c, d, a]);

areas = -atan2(Ca, spat(1)) - atan2(Cb, spat(2)) -atan2(Cc, spat( 3)) -atan2(Cd, spat(4));
% areas = mod(areas + 8.5 * pi, pi) - 0.5 * pi;

areas_func = matlabFunction(areas);
sines_func = matlabFunction(spat);


J_sines = jacobian(spat, reshape(x, 1, []));
J_area = reshape(gradient(areas, reshape(x, 1, [])), 3, 4);
J_area_sph = J_area - repmat(sum(J_area .* x), 3, 1) .* x;
J_sine_func_cell = {};

J_area_func = matlabFunction(J_area_sph);



load('cmx.mat');

areas_func(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4))
sines_func(cmx(1), cmx(2), cmx(3), cmx(4), cmx(5), cmx(6), cmx(7), cmx(8), cmx(9), cmx(10), cmx(11), cmx(12))
J_area = J_area_func(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4))

a1 = areas_func(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4));
cmx1 = cmx;
cmx1(2) = cmx1(2) + 1e-8;
cmx1(:, 1) = cmx1(:, 1) ./ norm(cmx1(:, 1));
a2 =  areas_func(cmx1(1, 1), cmx1(1, 2), cmx1(1, 3), cmx1(1, 4), cmx1(2, 1), cmx1(2, 2), cmx1(2, 3), cmx1(2, 4), cmx1(3, 1), cmx1(3, 2), cmx1(3, 3), cmx1(3, 4));
(a2 - a1) / 1e-8

% test 
data = load('workspace_1.mat');
m_x = data.m_x;
faces = data.faces;
clear data;
face_id = 6285;
cmx = m_x(faces(6285, :), :)';
J_area_sph = calc_area_jacobian_sphere(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4));
areas_func(cmx(1, 1), cmx(1, 2), cmx(1, 3), cmx(1, 4), cmx(2, 1), cmx(2, 2), cmx(2, 3), cmx(2, 4), cmx(3, 1), cmx(3, 2), cmx(3, 3), cmx(3, 4))

cmx1 = cmx;
cmx1(3) = cmx1(3) + 1e-9;
cmx1(:, 1) = cmx1(:, 1) ./ norm(cmx1(:, 1));
a2 =  areas_func(cmx1(1, 1), cmx1(1, 2), cmx1(1, 3), cmx1(1, 4), cmx1(2, 1), cmx1(2, 2), cmx1(2, 3), cmx1(2, 4), cmx1(3, 1), cmx1(3, 2), cmx1(3, 3), cmx1(3, 4));
(a2 - a1) / 1e-9

% output function to file 
% matlabFunction(J_area_sph, 'File', 'calc_area_jacobian_sphere.m', 'Comments', 'The function compute the area jacobian for a given face in the sphere');
matlabFunction(J_area_sph, 'File', 'calc_area_jacobian_sphere.m');

for i = 1 : 4
    J_sine_i = reshape(J_sines(i, :), 3, 4);
    J_sine_sph_i = J_sine_i - repmat(sum(J_sine_i .* x), 3, 1) .* x;
    function_name = sprintf('calc_sine_%d_jacobian_sphere', i);
    matlabFunction(J_sine_sph_i, 'File', function_name);
end

% save('area_sines_jacobian_functional_handle.mat', 'J_sine_func_cell', 'J_area_func');


