function [mesh_out] = reconstruct_spharm_descriptor_to_mesh(descriptor, components)
% the function aims to convert spharm descriptors to image. 
% it works for the cases: both cell and nuclear shape, only cell shape,
% only nuclear

max_deg = 31;
[vs, fs]=SpiralSampleSphere(4002);
% [vs fs] = sphereMesh([0 0 0 1]);
Zs = calculate_SPHARM_basis(vs, max_deg);

cell_descriptor = [];
nuc_descriptor = [];
if any(contains(components, 'cell'))
    cell_descriptor = descriptor(:, :, 1);
    if any(contains(components, 'nuc'))
        nuc_descriptor = descriptor(:, :, 2);
    end
else
    nuc_descriptor = descriptor(:, :, 1);
end

if ~isempty(cell_descriptor)
    fvec_cell = cell_descriptor;
    Zvert_cell = real(Zs * fvec_cell);
    mesh_out.cell_vertices = Zvert_cell;
    mesh_out.cell_faces = fs;
end

if ~isempty(nuc_descriptor)
    fvec_nuc = nuc_descriptor;
    Zvert_nuc = real(Zs * fvec_nuc);
    [Zvert_nuc, fs_nuc] = reorder_mesh_vertices(Zvert_nuc, fs, 'x-axis');
    mesh_out.nuc_vertices = Zvert_nuc;
    mesh_out.nuc_faces = fs_nuc;
end


end



