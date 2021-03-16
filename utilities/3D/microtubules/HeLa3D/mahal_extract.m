function mahal_dist = mahal_extract(feat_vector_MI,feat_vector_ori,in_cov,cellnum)

mahal_dist = (feat_vector_MI-feat_vector_ori) * (in_cov) * ...
        ((feat_vector_MI-feat_vector_ori)');
end
