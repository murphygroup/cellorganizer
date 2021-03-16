function synimage_run(i,k)
% Synthesize protein images using shape i and pattern k. The selected shape
% indices are [2,4,6,7,17]; the synthesized image is saved in ./newsamples/
% with name syn_i_k.

patternlist = {'LAM','Mit','Nuc','TfR'};

load(['./inter_results/image_synshapes/shape' int2str(k) '.mat'])
load(['./inter_results/protein_objects_gaussian/' ...
    patternlist{k} '_model'])

% An typical value of number of object is chosen N = [184 541 54 739];

%blank = nucimg(:,:,1);
%nucimg = cat(3,blank,nucimg,blank);

N = [184 541 54 739];

protimg = tp_genprotimage(nucimg,cellimg,model,N(i));
save(['./inter_results/image_synshapes/syn_' ...
    int2str(k) '_' int2str(i) '.mat'],'protimg')
