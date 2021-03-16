function [feat_vector_ori,idxes] = img2microtubule_feats(im_mt,cellmask, dnamask,centrosome_coordinate,w)
%feature finding function for microtubule model
%
%w is a binary vector of the length(methods) that indicates which features
%to compute
%
%merged the copies of getfeatvector in HeLa3D and HeLa2D 
%
%grj 10/4/14

% Get any images feature vector


if ~exist('centrosome_coordinate', 'var') | isempty(centrosome_coordinate)
    centrosome_coordinate = img2centrosome_coord(im_mt);
end

methods = {'haralick', ...
            'histpropwithcent', ...
            'harkd2', ...
            'harkd4', ...
            'totint', ...
            'radIntensity', ...
            'edge', ...
            'haralick2', ...
            'histpropwithcent2', ...
            'histbinwithcent2', ...
            'histDistLevel2', ...
            'statDistLevel', ...
            'distHist'};

if ~exist('w', 'var') | isempty(w)        
    w = ones(1,length(methods));
end
        
if length(w) ~= length(methods)
    error(['Variable ''w'' must be the same length as the number of methods (' num2str(length(methods)) ').']);
end
        
feat_vector_ori = [];

for i = 1:length(w)
    if w(i)
        idxes{i}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract2(im_mt,cellmask,dnamask, methods{i},centrosome_coordinate)];
        
        if any(isnan(feat_vector_ori))
            1;
        end
        
        idxes{i}(2) = size(feat_vector_ori,2);
    end
end


end % function ends
