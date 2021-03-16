function mainScript10_hela_acute(n_inputArray,rand_start,batchno,cellnum,subfolder)

% Loading Dboth from cell1 that has been semi automatically segmented
%load(['proc/cell_' num2str(cellnum) '/final_info.mat'],'Dbothfin','imgcent_coordinate','XYZres');
%load(['/home/jieyuel/lib_backup/Microtubules/proc/cell_' num2str(cellnum) '/final_info.mat'],'Dbothfin','imgcent_coordinate','XYZres');
load(['/proc/cell_' num2str(cellnum) '/final_info.mat'],'Dbothfin','imgcent_coordinate','XYZres');

folderName = ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n_inputArray)];
if ~exist(folderName,'dir')
mkdir(folderName);
end

% n_inputArray = [5,50:50:200];
%mulen_inputArray  = [5 10 15 20 25 30 35];
mulen_inputArray  = [5 10 15 20 25 30 35 40];
colli_inputArray = [0.9 0.95 0.98];

if subfolder==1  %%
coeff_var_Array = [0 0.1 0.2 0.3];
else  
coeff_var_Array = [0];  %%
end

rand_seed = rand_start;

tmpfolder = './tmp/';
if ~exist(tmpfolder,'dir')
   mkdir(tmpfolder);
end

for n = n_inputArray
  for mu_len = mulen_inputArray
    for coeffvar = coeff_var_Array
      sigma_len = coeffvar*mu_len;	
      for colli_min_number = colli_inputArray

       	  myfile = ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat'];
                tmpfile = [tmpfolder,regexprep(myfile,'/','_'),'.txt'];
                if ~exist(tmpfile,'file')
                   fid = fopen(tmpfile, 'w');
                   if (~exist(myfile,'file')) && (~exist([myfile(1:end-4),'-tcheck.mat'],'file'))
			if subfolder<=2
                      [K,tcheck] = colli_generator_acute(n,mu_len,sigma_len,colli_min_number,Dbothfin,imgcent_coordinate,XYZres,0.4,0.2,rand_seed,cellnum,subfolder,batchno);
            else
                rand_seed
			 [K,tcheck] = colli_generator_acute_jl2(n,mu_len,sigma_len,colli_min_number,Dbothfin,imgcent_coordinate,XYZres,0.4,0.2,rand_seed,cellnum,subfolder,batchno);
			end
                   end  % end of myfile
		     try
                   fclose(fid);
                   delete(tmpfile);
		     catch
		     end
                end  % end of tmpfile
	  
	  rand_seed = rand_seed + 1;

       end
     end
  end
end

 
end
