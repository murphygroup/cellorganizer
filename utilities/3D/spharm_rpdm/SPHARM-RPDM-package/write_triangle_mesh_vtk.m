function [fname] = write_triangle_mesh_vtk(faces, vertices, fname)
% write triangle mesh to vtk format 
% xruan: 03/11/2018

% write header 

nvert = size(vertices, 1);
nface = size(faces, 1);

hd_string = sprintf('# vtk DataFile Version 4.1\nvtk output\nASCII\nDATASET POLYDATA\n');
hd_string = sprintf('%sPOINTS %d float\n', hd_string, nvert);

fid = fopen(fname, 'w');
fwrite(fid, hd_string);
fclose(fid);

vertices = vertices';
term_ind = floor(nvert / 3) * 3;
vertices_1 = reshape(vertices(:, 1 : term_ind), 9, [])';
dlmwrite(fname, vertices_1, '-append', 'delimiter', ' ');
if term_ind < nvert
    vertices_2 = reshape(vertices(:, term_ind + 1 :end), 1, []);
    dlmwrite(fname, vertices_2, '-append', 'delimiter', ' ');
end

polygon_info = sprintf('POLYGONS %d %d\n', nface, nface * 4);
fid = fopen(fname, 'a');
fwrite(fid, polygon_info);
fclose(fid);

faces = [ones(nface, 1) * 3, faces-1];
dlmwrite(fname, faces, '-append', 'delimiter', ' ');

end















