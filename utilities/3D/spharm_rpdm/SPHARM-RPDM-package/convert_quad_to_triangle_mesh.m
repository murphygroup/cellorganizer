function [ faces ] = convert_quad_to_triangle_mesh( faces )
% convert quad mesh to triangle mesh. 
% based on the code from SPHARM-MAT
% xruan: 03/11/2018

if size(faces,2)==4
    qfaces = faces;
    faces = [faces(:,1:3); faces(:,[3 4 1])];
    
    dif1 = faces(:,1)-faces(:,2); dif2=faces(:,2)-faces(:,3); dif3=faces(:,3)-faces(:,1);
    indDif1 = find(dif1 == 0);
    indDif2 = find(dif2 == 0);
    indDif3 = find(dif3 == 0);
    % indDif = union(indDif1, indDif2, indDif3);
    % 01/22/2018 xruan
    indDif = unique([indDif1; indDif2; indDif3]);
    ufIDX = setdiff([1:size(faces,1)], indDif);
    faces = faces(ufIDX,:);
end


end

