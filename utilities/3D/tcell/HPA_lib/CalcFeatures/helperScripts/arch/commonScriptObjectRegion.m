for tmp_i=1:length(names)
    names_l{tmp_i} = [names_l{tmp_i} '_large'];
end

feats = [feats feats_l];
names = [names names_l];
slfnames = [slfnames slfnames];
