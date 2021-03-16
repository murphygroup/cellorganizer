function [processed_vertex_values] = process_mesh_neighborhoods(single_mesh, vertex_values, compute_function)
  % Apply operations on the discrete neighborhood of each vertex (like convolution or morphological image processing, but on meshes).
  %
  % 2013-02-27 tebuck: created.
  %
  % Dependencies:
  % From the File Exchange: toolbox_graph
  adjacency = logical(triangulation2adjacency(single_mesh.faces, single_mesh.vertices));
  processed_vertex_values = zeros(size(adjacency, 1), 1);
  for vertex_index = 1:length(processed_vertex_values)
    processed_vertex_values(vertex_index) = compute_function(vertex_values(vertex_index), vertex_values(adjacency(vertex_index, :)));
  end
end