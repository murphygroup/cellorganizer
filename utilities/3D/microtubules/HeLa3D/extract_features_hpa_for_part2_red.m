function extract_features_hpa_for_part2_red(subfolder,cellnum,batchno,n_inputArray,mulen_inputArray,colli_inputArray,coeff_var_Array)

setgenpath
addpath slicfiles/
imCategory = 'general';

featTypeList{1} = 'haralick';
featTypeList{2} = 'histpropwithcent';
featTypeList{3} = 'harkd2';
featTypeList{4} = 'harkd4';
featTypeList{5} = 'totint';
featTypeList{6} = 'radIntensity';
featTypeList{7} = 'edge';

  featTypeList{8} = 'haralick2';
  featTypeList{9} = 'histpropwithcent2';
  featTypeList{10} = 'histbinwithcent2';
  featTypeList{11} = 'histDistLevel2';
  featTypeList{12} = 'statDistLevel';
  featTypeList{13} = 'distHist';

%subfolder = 1;

featType = 'all';

fin_mat = [];
I = 1;

tmpfolder = './tmp/';
if ~exist(tmpfolder,'dir')
   mkdir(tmpfolder);
end

savefile = ['outputs_' num2str(subfolder) '/featvals/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/' imCategory '/' featType '/forpart2_all_feats_G_psf_batch_' num2str(batchno) '.mat'];
tmpfile = [tmpfolder,regexprep(savefile,'/','_'),'.txt'];
if ~exist(tmpfile,'file')
   fid = fopen(tmpfile,'w');
   if ~exist('savefile','file')

[protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(cellnum);

libpars = reducedlibrary(cellnum,subfolder);

for I = 1:size(libpars,1)
        n = libpars(I,1);
        mu_len = libpars(I,2);
        sigma_len = libpars(I,3);
        colli_min_number = libpars(I,4);

	myfile = ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat'];

	if (exist(myfile,'file') == 2)
		[G_psf] = getsynimage_hela(n,mu_len,sigma_len,colli_min_number,cellnum,batchno,subfolder);
	 
		feat_vector = [];
		%for J = 1:7
		for J = 1:13
			feat_vector = [feat_vector,feat_extract(G_psf,[],featTypeList{J},imgcent_coordinate,cellnum)];
		end
		fin_mat(I,:) = [n mu_len sigma_len colli_min_number feat_vector];
		%I = I + 1 
	end
end

if ~exist(['outputs_' num2str(subfolder) '/featvals/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/' imCategory '/' featType '/'],'dir')
mkdir(['outputs_' num2str(subfolder) '/featvals/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/' imCategory '/' featType '/']);
end

save(['outputs_' num2str(subfolder) '/featvals/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/' imCategory '/' featType '/forpart2_all_feats_G_psf_batch_' num2str(batchno) '.mat'],'fin_mat');

end % end of savefile
try
fclose(fid);
delete(tmpfile);
catch
end
end % end of tmpfile

end % end of function
