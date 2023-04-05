function models2report_v2(models, param, classlabels, fileID)
% MODELS2_V2 Generate a report comparing various SLML models.
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% models                   A cell array of SLML models
% param                    A structure holding the function options
% classlabels
% fileID

% List Of Parameters        Descriptions
% ------------------        ------------
% targetDirectory           (optional) Directory where the images are going to be saved. Default is current directory.
% prefix                    (optional) Filename prefix for the synthesized images. Default is 'demo'
% numberOfSynthesizedImages (optional) Number of synthesized images. Default is 1.
% compression               (optional) Compression of tiff, i.e. 'none', 'lzw' and 'packbits'
% debug                     (optional) Keeps temporary results and catches errors with full reports
% display                   (optional) Display flag for figures
% verbose                   (optional) Print the intermediate steps to screen. Default is true.
% microscope                (optional) Microscope model from which we select a point spread function. Default is 'none'
% synthesis                 (optional) Synthesis parameter that allows to
% synthesize                           'nucleus', 'framework' or 'all'. Default is 'all'
% sampling.method           (optional) Can be 'disc' or 'sampled'. Default is trimmed
% sampling.method.density   (optional) An integer. Default is empty.
% protein.cytonuclearflag   (optional) Can 'cyto', 'nucleus' or 'all'. Default is all.
% resolution.cell           (optional) The resolution of the cell and nucleus that are being passed in
% resolution.objects        (optional) The resolution of the object model being synthesized
% image_size                (optional) The image size. Default is [1024 1024] for both 2D and 3D in x and y
%
% output.tifimages           (optional) boolean flag specifying whether to write out tif images
% output.indexedimage        (optional) boolean flag specifying whether to write out indexed image
% output.blenderfile         (optional) boolean flag specifying whether to write out (.obj) files for use in blender
% output.blender.downsample  (optional) downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size)

%
% Copyright (C) 2018 Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu
%
% Feb. 16, 2023 R.F.Murphy correct test for whether nuclear shape model is
% SPHARM; save comparison scores to .mat file
% March 21, 2023 R.F.Murphy add shape space plot to spharm-obj


%put includespharm under the includeprot
param = ml_initparam(param, ...
    struct('verbose',true,'includenuclear',true, ...
    'includecell',true,'includeprot',true));

n = length(models);

indvarlabels = classlabels;
legendentry = regexp(indvarlabels, ';', 'split');

colors = jet(n);

is3D = strcmpi(models{1}.dimensionality, '3D');

[ nuclearparams, cellparams, proteinparams] = getModelReportParams( models{1}, param );

protfield = 'proteinModel';
indvar = 1:n;

%%%%%%%%%%%%%%%%%%%%%% MODEL IDENTIFIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html(fileID, 'Model Identifiers');
for i=1:n
    identifiers = {'filename' 'name' 'id'};
    sout = [];
    for j=1:3
        identifier = identifiers{j};
        model = models{i};
        if isfield(model, identifier)
            identifier_value = eval(sprintf('model.%s', identifier));
            sout = strcat(sout,sprintf( 'model%d.%s = %s;', i, identifier, identifier_value ));
        else
            sout = strcat(sout,sprintf( 'model%d.%s = %s', i, identifier, 'UNSET'));
        end
    end
    text2html(fileID,sout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compare all model parameters
% first collect protein parameters
y = []; 
if param.includeprot 
    for j = 1:n
      if is3D && strcmpi(models{1}.proteinModel.type,'gmm')
         compartments{j} = eval(['models{j}.' protfield '.cytonuclearflag']);
      else
        compartments{j} = 'all';
      end
    end
     
    for j=1:n
        for i=1:length(proteinparams)
            y(j, i)=eval(['models{j}.' protfield '.' proteinparams{i}]');
        end
    
%check if gmm
        if strcmpi(models{1}.proteinModel.type,'gmm')       
            if is3D
                [positionBeta, fractdist] = beta2posMap(models{j}.proteinModel.position.beta);
            else
                [positionBeta, fractdist] = beta2posMap(models{j}.proteinModel.positionModel.beta);
            end
        
            objpdf = squeeze(sum(sum(positionBeta,1),3));

            if strcmpi(compartments{j}, 'cyto')
                domain = fractdist >= 0;
                objpdf(~domain) = [];
                objpdf = [zeros(1, 1+sum(~domain)), objpdf ] ;
                fractdists{j} = [fractdist(~domain) 0 fractdist(domain)];
        
            elseif strcmpi(compartments{j}, 'nuc')
                domain = fractdist <= 0;
                objpdf(~domain) = [];
                objpdf = [objpdf zeros(1, 1+sum(~domain))];
                fractdists{j} = [fractdist(domain) 0 fractdist(~domain)];
            
            else
                fractdists{j} = fractdist;
            end
        
            %renormalize
            objpdf = objpdf./sum(objpdf);
        
            objpdfs{j} = objpdf;
        end
     end 
  end
%
% now compare nuclear parameters
%
 y2 = []; 
 y3 = []; cellparams = [];
 if param.includenuclear   
       
        for i=1:length(nuclearparams)
            for j=1:n
                y2(j,i)=eval(['models{j}.nuclearShapeModel.' nuclearparams{i}]');
            end
        end

 end

 
 paramnames = [proteinparams nuclearparams cellparams];
 yall = [y y2 y3];

%
% plot all of the parameters for all of the models individually, and plot
% the ones that change the most together
%

%%%%%%%%%%%%%%%%%%%%%% PARAMETER COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(yall)>0 
  model1_name = models{1}.name;
  model2_name = models{2}.name;
  fig_count = 0;
  %fig_names = {sprintf('Protein Model Parameters'), ...
  %  sprintf('Nucleus and Cell Model Parameters'), ...
  %  'Coefficient of Variation', ...
  %  'Fractional Distance', ...
  %  'Number of Objects'};
  figures = [];
  npanel = length(paramnames);
  ncolumns = 4; nrows = 3; ntot = ncolumns*nrows;
  for iplot=1:ceil(npanel/ntot)
    figures = [figures, figure];
    fig_count = fig_count + 1;
    nleft = npanel-(iplot-1)*ntot;
    for i=(iplot-1)*ntot+1:(iplot-1)*ntot+min(nleft,ntot)
        subplot((ceil(min(nleft,ntot)/ncolumns)),ncolumns,i-(iplot-1)*ntot);
        plot(indvar,yall(:,i),'rd-');
        xlabel(indvarlabels);
        ylabel(paramnames{i});
    end
    saveas(figures(fig_count), sprintf('image%d.png', fig_count), 'png');
    fP = figures(fig_count).Position;
    figures(fig_count).Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    saveas( figures(fig_count), sprintf('image%d_thumbnail.png', fig_count), 'png' ); 3/21/2023
    img2html(fileID, sprintf('image%d.png', fig_count), ...
        sprintf('image%d_thumbnail.png', fig_count), 'Parameter comparison');
  end

%%%%%%%%%%%%%%%%%%%%% PARAMETERS WITH LARGEST CHANGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 figures = [figures, figure];
 fig_count = fig_count + 1;
 normy = yall./repmat(mean(yall),n,1);
 cv=std(normy);
 [sorted,order]=sort(cv,'descend');
 % find the ones whos CV is greater than 50%
 idx = find(sorted>0.5);
 mintodisplay = min(3,length(order));
 if length(idx) > mintodisplay
    order = order(sorted>0.5);
 else
    order = order(1:mintodisplay);
 end

 colorstr = {'rd-','g+-','b^-','co-','mv-','yx-'};
 for i=1:min(length(colorstr),length(order))
    plot(indvar,normy(:,order(i)),colorstr{i});
    legends{i} = paramnames{order(i)};
    hold on;
 end
 legend(legends);
 xlabel(indvarlabels);
 ylabel('Coefficient of variation');
 title('Parameters ordered by extent of variation');
 hold off;

 saveas(figures(fig_count), sprintf('image%d.png', fig_count), 'png');
 fP = figures(fig_count).Position;
 figures(fig_count).Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
 saveas( figures(fig_count), sprintf('image%d_thumbnail.png', fig_count), 'png' ); 3/21/2023
 img2html(fileID, sprintf('image%d.png', fig_count), ...
     sprintf('image%d_thumbnail.png', fig_count), 'Parameters with largest differences');
end

%%%%%%%%%%%%%%%%%%%%% GMM PROTEIN MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.includeprot
    if strcmpi(models{1}.proteinModel.type,'gmm')
        figures = [figures, figure];
        fig_count = fig_count + 1;
        for i = 1:length(models)
        %marginalize out the anglular dimensions.

            plot(fractdists{i}, objpdfs{i} , 'color', colors(i,:), 'LineWidth',2)
            hold on;
        end

        xlabel('Fractional distance between nuclear and plasma membranes')
        ylabel('Relative probability density')
        legend(legendentry);
        set(gca,...
        'YTickLabel','')

        saveas(figures(fig_count), sprintf('image%d.png', fig_count), 'png');
        fP = figures(fig_count).Position;
        figures(fig_count).Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
        saveas( figures(fig_count), sprintf('image%d_thumbnail.png', fig_count), 'png' ); 3/21/2023
        img2html(fileID, sprintf('image%d.png', fig_count), ...
            sprintf('image%d_thumbnail.png', fig_count), 'GMM Vesicle Model Comparison');

%%%%%%%%%%%%%%%%%%%%% PROTEIN GMM SPATIAL DISTRIBUTION COMPARISON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if is3D
        

            objmu = y(:,strcmpi(proteinparams, 'frequency.mu'));
            objsigma = y(:,strcmpi(proteinparams, 'frequency.sigma'));
            xrange = [0 max(exp(objmu+objsigma))];

            figures=[figures, figure];
            for i = 1:length(models)
                semilogx(xrange(1):xrange(2), lognpdf(xrange(1):xrange(2), objmu(i), objsigma(i)), 'color', colors(i,:), 'LineWidth',2)
                hold on
            end
            xlabel('Number of objects')
            ylabel('Relative probability density')
            legend(legendentry);
            set(gca,...
            'YTickLabel','')
        else
            objalpha = y(:,strcmpi(proteinparams, 'objectModel.numStatModel.alpha'));
            objbeta = y(:,strcmpi(proteinparams, 'objectModel.numStatModel.beta'));

            xrange = 0:1000;
        
            for i = 1:length(models)
                gamcdfs(i,:) = gamcdf(xrange, objalpha(i), objbeta(i));
                gampdfs(i,:) = gampdf(xrange, objalpha(i), objbeta(i));
            end
        
            inds = any(gamcdfs < 0.95,1);
        
            figures=[figures, figure];
        
            for i = 1:length(models)
                plot(xrange(inds), gampdfs(i,inds)','color', colors(i,:))
                legendentry{i} = num2str(indvar(i));
                hold on
            end
        
            xlabel('Number of objects')
            ylabel('Relative probability density')        
            legend(legendentry);
            set(gca,...
            'YTickLabel','')
        end
        fig_count = fig_count + 1;
        saveas(figures(fig_count), sprintf('image%d.png', fig_count), 'png');
        fP = figures(fig_count).Position;
        figures(fig_count).Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
        saveas( figures(fig_count), sprintf('image%d_thumbnail.png', fig_count), 'png' ); 3/21/2023
        img2html(fileID, sprintf('image%d.png', fig_count),...
            sprintf('image%d_thumbnail.png',fig_count), 'Spatial model comparison');
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPHARM_OBJ_MODEL COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(models{1}.proteinModel.type,'spharm_obj')
        header2html(fileID, 'spharm-obj Model Comparison');

        for i=1:length(models)
            models_ShapeModel{i} = models{i}.proteinModel.spharm_obj_model.cellShapeModel;
        end
%        [consolidated_model,shapecomp,shapepvals]=compare_shape_space_spharm_obj(models_ShapeModel,fileID,param);
%        total_variation=consolidated_model.total_variation;
%        text2html(fileID, ['Total_Variation=''%s'';\n', total_variation]);

        % calculate shape comparison statistics 3/21/2023
        [reducedSPHARMVec,shapecomparisonscores,shapecomparisonpvals]=compare_shape_model_spharm_obj(models_ShapeModel, param, fileID);

        % generate composite shape space figure 3/21/2023
        header2html(fileID, 'spharm-obj Organelle Shape Model Comparison');
        text2html(fileID, sprintf('Number of objects: Model 1 %i; Model 2 %i',...
            size(models_ShapeModel{1}.cell_params,3),...
            size(models_ShapeModel{2}.cell_params,3)));
        compare_shape_space_spharm_obj(models_ShapeModel,fileID,param);
        
        for j=1:length(models)-1
        %spatial model comparision
            models_SpatialModel{1}=models{j}.proteinModel.spharm_obj_model.spatial;
            models_SpatialModel{2}=models{end}.proteinModel.spharm_obj_model.spatial;

            header2html(fileID, 'Spatial Model Comparison');
    %       compare_spatial_model_spharm_obj(models_SpatialModel,fileID,param);
            %clique percolation
            d1 = length(models_SpatialModel{1}.normdists);
            d2 = length(models_SpatialModel{2}.normdists);

            normdists=horzcat(models_SpatialModel{1}.normdists, models_SpatialModel{2}.normdists);
            anglestheta=horzcat(models_SpatialModel{1}.anglestheta, models_SpatialModel{2}.anglestheta);
            anglesphi=horzcat(models_SpatialModel{1}.anglesphi, models_SpatialModel{2}.anglesphi);

            [x,y,z]=sph2cart(anglestheta,anglesphi,normdists);
            
            pos = horzcat(x',y',z');
            [coeff,score,latent] = pca(pos);
            
            prob2 = d1/(d1+d2);
            temp = zeros(1, d1+d2);
            for i=1:(d1+d2)
                Idx = knnsearch(pos, pos(i, :), 'K', 9);
                prob1 = sum(Idx(2:end) < d1)/8;
                temp(i)=prob1.*log((prob1 + 1e-7)/prob2) + (1-prob1).*log((1-prob1+ 1e-7)/(1-prob2));
            end

            spatialcomparisonscores(j) = mean(temp);
            text2html(fileID, sprintf('Spatial distribution comparison score between model%d and model%d: %.4f;\n', j, length(models), mean(temp)));
            
            KL_div = show_spatial_distribution(models_SpatialModel,fileID);
            p_val = spatial_permutation_test(KL_div, x(d1+1:end), y(d1+1:end), z(d1+1:end));
            spatialcomparisonspval(j) = p_val;
            text2html(fileID, sprintf('P value of spatial distribution divergence: %.4f;\n', p_val));

            %p_val = permutation_test(mean(temp), pos(d1+1:end,:));
        end

        X = [ones(length(reducedSPHARMVec{end}), 1) pos(length(reducedSPHARMVec{j})+1:end, :)];
        B = X\reducedSPHARMVec{end}(:, 1);
        yCalc = X*B;
        shapecomparisonsR2 = 1 - sum((reducedSPHARMVec{end}(:, 1) - yCalc).^2)/sum((reducedSPHARMVec{end}(:, 1) - mean(reducedSPHARMVec{end}(:, 1))).^2);
        text2html(fileID, sprintf('R2: %.4f;\n', shapecomparisonsR2));

        save('spharm_obj_comparison','shapecomparisonscores','shapecomparisonpvals',...
            'spatialcomparisonscores','spatialcomparisonspval','shapecomparisonsR2')
    end
end

if param.includecell
    if strcmpi(models{1}.cellShapeModel.type,'spharm_rpdm')
        models_cellShapeModel{1} = models{1}.cellShapeModel;
        models_cellShapeModel{2} = models{2}.cellShapeModel;
        header2html(fileID, 'SPHARM-RPDM Cell Shape Model Comparison');
        text2html(fileID, sprintf('Number of cells: Model 1 %i; Model 2 %i',...
            size(models{1}.cellShapeModel.cell_params,3),...
            size(models{2}.cellShapeModel.cell_params,3)));
        compare_shape_space_spharm_obj(models_cellShapeModel,fileID,param);
    end
end
 
if param.includenuclear
    if strcmpi(models{1}.nuclearShapeModel.type,'spharm_rpdm') %corrected 2/13/2023
        models_nucShapeModel{1} = models{1}.nuclearShapeModel;
        models_nucShapeModel{2} = models{2}.nuclearShapeModel;
        compare_shape_space_spharm_obj(models_nucShapeModel,fileID,param);
    end
end
 
end

