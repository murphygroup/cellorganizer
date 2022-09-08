function [mesh_out] = reconstruct_spharm_descriptor_to_mesh(descriptor, components)
% the function aims to convert spharm descriptors to image. 
% it works for the cases: both cell and nuclear shape, only cell shape,
% only nuclear

% 8/20/2022 Ted make this function deployable as a binary to be used in
% Docker containers

% if isdeployed
%     disp('Running deployed version of reconstruct_spharm_descriptor_to_mesh...');
% 
%     %getting info read into matlab
%     %when method is deployed
%     if length(varargin) == 1
%         text_file = varargin{1};
%         
%         [filepath, name, ext] = fileparts(text_file);
% 
%         if ~exist(text_file, 'file')
%             warning('Input file does not exist. Exiting method.');
%             return
%         end
% 
%         disp(['Attempting to read input file ' text_file]);
%         fid = fopen(text_file, 'r' );
% 
%         disp('Evaluating lines from input file');
%         while ~feof(fid)
%             line = fgets(fid);
%             disp(line);
%             try
%                 eval(line);
%             catch err
%                 disp('Unable to parse line');
%                 getReport(err)
%                 return
%             end
%         end
%         fclose(fid);
%         if ~exist('options', 'var')
%             options = {};
%         end
% 
%         model = load(model_path);
%         descriptor = model.cellShapeModel.all_spharm_descriptors;
%         components = model.cellShapeModel.components;
% 
%         cellind = options.cellind;
% 
%         if numel(model.components) > 1
%             descriptors = reshape(descriptors, [], 2, size(descriptors, 2), size(descriptors, 3));
%             descriptors = permute(descriptors, [1, 3, 2, 4]);
%             descriptors = descriptors(: , :, :, cellind);
%         end
% 
%         descriptors = descriptors(: , :, cellind);
%  
%     else
%         descriptor = varargin{1};
%         components = varargin{2};
%     end
%     
%     
% else
%     
%     descriptor = varargin{1};
%     components = varargin{2};
%     
% end

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

%save mesh if deployed
% if isdeployed
%     disp('saving mesh...');
%     output_dir = join([options.output_dir, '/mesh_output.mat']);
%     save(output_dir, 'mesh_out');
% 
% end

end



