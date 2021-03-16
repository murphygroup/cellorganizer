function [gausssize,gaussinten,objnum, options]= gmm_objempdistr_percell(mixes, objintens,options)
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
%
% June 13, 2013 D. Sullivan  Refactored code to compute model parameters in
%                            a per-cell manner for parallel computing.

%protype = {'LAM','Nuc','Mit','TfR'};
% files = ml_dir([gaussobjpath filesep '*gaussobjs.mat']);

%for each pattern in the protype list
% for pidx = 1:1:length(files)

    % Intensity fit and number fit
    gaussinten = [];
    gausssize = [];

%     for i = 1:length(mixes)
        objnum = 0;
        for j = 1:length(mixes{1})
            %if ~isempty(mixes{i}(j))
               %if ~isempty(mixes{i}(j).ncentres)
                for k = 1:mixes{1}(j).ncentres
                    %disp( [ '(i,j,k)=(' num2str(i) ',' num2str(j) ',' num2str(k) ')']);
                    gaussinten = [gaussinten;...
                        objintens{1}(j)*mixes{1}(j).priors(k)];
                    sigma = mixes{1}(j).covars(:,:,k);
                    
                    %Z dimension - resize z to make cubic voxels. This is
                    %necessary 
                    if isfield(options.model,'proteinUpsampleZ')
                        sigma(:,end) = sigma(:,end) * options.model.proteinUpsampleZ;
                        sigma(end,:) = sigma(end,:) * options.model.proteinUpsampleZ;

                    end
                    
                    lambda = sqrt(eig(sigma));
                    
                    if ~all(isreal(lambda))
                        warning('Removing imaginary GMM eigen-value components.')
                        lambda = real(lambda);
                    end
                    
                    
                    gausssize = [gausssize;fliplr(lambda')];
                   %catch
                    %disp( ['Ignoring object ' num2str(k) '.' ] );
                   %end
                end
               %end
               %end                
%                 try
                 objnum = objnum + mixes{1}(j).ncentres;
%                 catch
                 
%                 end
        end
%         end    
    
     
% end

%D. Sullivan 2/24/13 Adjust the protein_resolution given the upsamplings
% %X dimension
% if isfield(param,'proteinUpsampleX')
%     param.model.protein_resolution(1) = param.model.protein_resolution(1)*param.proteinUpsampleX;
% end
% %Y dimension
% if isfield(param,'proteinUpsampleY')
%     param.model.protein_resolution(2) = param.model.protein_resolution(2)*param.proteinUpsampleY;
% end
if isfield(options.model,'proteinUpsampleZ')
    options.model.protein_resolution(3) = options.model.protein_resolution(3)/options.model.proteinUpsampleZ;
end
