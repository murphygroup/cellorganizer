function nvardist = nvardist_extract(feat_vector_MI,feat_vector_ori,varmat)

nvardist = sum(((feat_vector_MI-feat_vector_ori).^2)./varmat);
end
