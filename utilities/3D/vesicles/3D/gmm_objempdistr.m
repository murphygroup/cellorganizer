function param = gmm_objempdistr(gaussobjpath,savepath,param)
% Combining single object property for statistical analysis including
% object size, number, and intensities
%
%

%Author: 
%
%Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% Feb 22, 2013 D. Sullivan   Removed the fixed scaling of z by 5. 
% Feb 24, 2013 D. Sullivan   Added param struct to pass in the upsampling
%                            of the z dimension of the protein pattern.
%                            param.proteinUpsampleZ. Returns updated model
%                            resolution in param structure

%protype = {'LAM','Nuc','Mit','TfR'};
files = ml_dir([gaussobjpath filesep '*gaussobjs.mat']);

%for each pattern in the protype list
for pidx = 1:1:length(files)
  if ~exist([savepath filesep files{pidx}], 'file' )

    %load the matfile that contains the gaussian objects
    load([gaussobjpath filesep files{pidx}]);

    % Intensity fit and number fit
    intensities = [];
    objsize = [];

    for i = 1:length(mixes)
        objnum(i) = 0;
        for j = 1:length(mixes{i})
            %if ~isempty(mixes{i}(j))
               %if ~isempty(mixes{i}(j).ncentres)
                for k = 1:mixes{i}(j).ncentres
                    %disp( [ '(i,j,k)=(' num2str(i) ',' num2str(j) ',' num2str(k) ')']);
                    intensities = [intensities;...
                        objintens{i}(j)*mixes{i}(j).priors(k)];
                    sigma = mixes{i}(j).covars(:,:,k);
                    
                    %D. Sullivan 2/22/13 rather than blindly multiplying by
                    %5, see findobjs_run.m where the resolution of the
                    %objects is pre-adjusted s.t. each dimension has the
                    %same resolution
%                     sigma(:,end) = sigma(:,end) * 5;
%                     sigma(end,:) = sigma(end,:) * 5;
                    
                    %D. Sullivan 2/24/13 
                    %Reverting to post processing adjustments but using a
                    %param field now instead of just "5". Resolution is
                    %appropriately adjusted. for each resolution
                    %Note, X and Y dimensions are commented out for now
                    %unitl they are tested.
%                     %X dimension
%                     if isfield(param,'proteinUpsampleX')
%                         sigma(:,1) = sigma(:,end) * param.proteinUpsampleX;
%                         sigma(1,:) = sigma(end,:) * param.proteinUpsampleX;
%                     end
%                     %Y dimension
%                     if isfield(param,'proteinUpsampleY')
%                         sigma(:,2) = sigma(:,end) * param.proteinUpsampleY;
%                         sigma(2,:) = sigma(end,:) * param.proteinUpsampleY;
%                     end
                    %Z dimension
                    if isfield(param.model,'proteinUpsampleZ')
                        sigma(:,end) = sigma(:,end) * param.model.proteinUpsampleZ;
                        sigma(end,:) = sigma(end,:) * param.model.proteinUpsampleZ;

                    end


                    
                    lambda = sqrt(eig(sigma));
                    objsize = [objsize;fliplr(lambda')];
                   %catch
                    %disp( ['Ignoring object ' num2str(k) '.' ] );
                   %end
                end
               %end
               %end                
                try
                 objnum(i) = objnum(i) + mixes{i}(j).ncentres;
                catch
                 
                end
            end
        end    
    
     if ~exist( savepath, 'dir' )
       mkdir( savepath );
     end
      save([savepath filesep files{pidx}],'objsize','intensities','objnum')
 else  
disp( 'Intermediate results found, skipping recalculation' );
 end
end

%D. Sullivan 2/24/13 Adjust the protein_resolution given the upsamplings
% %X dimension
% if isfield(param,'proteinUpsampleX')
%     param.model.protein_resolution(1) = param.model.protein_resolution(1)*param.proteinUpsampleX;
% end
% %Y dimension
% if isfield(param,'proteinUpsampleY')
%     param.model.protein_resolution(2) = param.model.protein_resolution(2)*param.proteinUpsampleY;
% end
if isfield(param.model,'proteinUpsampleZ')
    param.model.protein_resolution(3) = param.model.protein_resolution(3)/param.model.proteinUpsampleZ;
end
