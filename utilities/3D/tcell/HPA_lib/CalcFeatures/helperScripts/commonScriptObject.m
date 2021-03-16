
feats = zeros(1,18);
if optimize && protein_channel_blank
  feats(:) = nan;
else
  
feats(1) = mean(protobjs);
feats(2) = var(protobjs);
try
feats(3) = max(protobjs)/min(protobjs);
catch
%keyboard
if isempty(protobjs)
   feats(3)=0;
else
   error('Check protobjs!');
end
end
feats(4) = length(protobjs)/length(nucobjs);
feats(5) = normEuN;

if isempty(protobjs_mthr)
    feats([6:10 15:18]) = nan;
else
    feats(6) = mean(protobjs_mthr);
    feats(7) = var(protobjs_mthr);
    feats(8) = max(protobjs_mthr)/min(protobjs_mthr);
    feats(9) = length(protobjs_mthr)/length(nucobjs);
    feats(10) = normEuN_mthr;

    feats(15) = mean(largeprotobjs_mthr);
    feats(16) = var(largeprotobjs_mthr);
    feats(17) = max(largeprotobjs_mthr)/min(largeprotobjs_mthr);
    feats(18) = length(largeprotobjs_mthr)/length(nucobjs);
end

feats(11) = mean(largeprotobjs);
feats(12) = var(largeprotobjs);
feats(13) = max(largeprotobjs)/min(largeprotobjs);
feats(14) = length(largeprotobjs)/length(nucobjs);

end

names = {'Avg. # pixels per obj.','Var. of # pixels per obj.','Ratio of size of largest obj. to smallest', 'nuc: # prot objs to nuc objs','normalized Euler number', ...
         'Avg. # pixels per obj._mthr','Var. of # pixels per obj._mthr','Ratio of size of largest obj. to smallest_mthr', 'nuc: # prot objs to nuc objs_mthr','normalized Euler number_mthr', ...
    'Avg. # pixels per large obj.','Var. of # pixels per large obj.','Ratio of size of largest large obj. to smallest', 'nuc: # large prot objs to nuc objs', ...
    'Avg. # pixels per large obj._mthr','Var. of # pixels per large obj._mthr','Ratio of size of largest large obj. to smallest_mthr', 'nuc: # large prot objs to nuc objs_mthr'};
slfnames = {'SLF1.3','SLF1.4','SLF1.5','SLF1.2','' ,'','','',''};
