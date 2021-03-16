function [result_image] = robust_read_stack(given_filename)
  % Given one or a cell array of possible images' filenames, find one that is a 3D stack of images and return it:
  % 
  % 2013-11-08 tebuck: Copied from master_script_produce_windows.m.
  % 
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  

  all_given_filenames = given_filename;
  if ~iscell(given_filename)
    all_given_filenames = {given_filename};
  end
  result_image = [];
  for given_filename_index = 1:length(all_given_filenames)
    given_filename = all_given_filenames{given_filename_index};
    if ~exist(given_filename, 'file')
      continue
    end
    [given_filename_directory, given_filename_base, given_filename_extension] = fileparts(given_filename);
    image_information = imfinfo(given_filename);
    number_slices = length(image_information);
    result_image = cell2mat(reshape(arrayfun(@(x)x.data, tiffread27(given_filename, [], false, false), 'UniformOutput', false), 1, 1, []));
    if size(result_image, 3) == 1
      result_image = zeros([image_information(1).Height, image_information(1).Width, number_slices]);
      for slice_index = 1:number_slices
        result_image(:, :, slice_index) = getfield(tiffread27(given_filename, slice_index, false, false), 'data');
      end
    end
    if size(result_image, 3) > 1
      break
    end
  end
  result_image = double(result_image);
  if isempty(result_image) || size(result_image, 3) <= 1
    error('None of the image files given exist, are non-empty, or contain more than one Z slice!')
  end
  
  
end
  
