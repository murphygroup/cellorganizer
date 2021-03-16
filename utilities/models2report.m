function figures=models2report( models, param, classlabels, fileID)
% MODELS2REPORT Generate a report comparing various SLML models.
%
% COMMENT: When the method is not deployed the shape of the input arguments is
% slml2report( models, param )
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% models                   A cell array of filenames
% param                    A structure holding the function options
%
% The shape of param is described
%
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
% Example
% -------
% instances = { 'model01.xml', 'model02.mat', 'model03.mat' };
% param.targetDirectory = pwd;
% param.prefix = 'demo'
% param.numberOfSynthesizedImages = 100;
% param.compression = 'lzw';
% param.microscope = 'svi';
% param.verbose = true;
%
% >> slml2img( instances, param );

% Author: Robert F. Murphy
% Created: June 29, 2013

% Copyright (C) 2013-2016 Murphy Lab
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

% if ~isdeployed
%     if nargin ~= 2
%         error( ['CellOrganizer: Wrong number of input arguments. If slml2img ' ...
%             ' is not deployed, then the number of input arguments must be 2.'] );
%     end
%
% end

%INPUT ARGUMENTS
param = ml_initparam(param, ...
    struct('verbose',true,'includenuclear',true, ...
    'includecell',true,'includeprot',true));

n = length(models);

indvarlabels = classlabels;
legendentry = regexp(indvarlabels, ';', 'split');

colors = jet(n);

is3D = strcmpi(models{1}.dimensionality, '3D');

[ nuclearparams, ~, proteinparams ] = getModelReportParams( models{1}, param );

protfield = 'proteinModel';
indvar = 1:n;

%
% compare protein model parameters
%
for j = 1:n
    if is3D
    compartments{j} = eval(['models{j}.' protfield '.cytonuclearflag']);
    else
        compartments{j} = 'all';
    end
end

y = []; 
if param.includeprot  
    for j=1:n
        for i=1:length(proteinparams)
            y(j, i)=eval(['models{j}.' protfield '.' proteinparams{i}]');
        end
        
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
%
% compare nuclear parameters
%
y2 = []; 
y3 = []; cellparams = [];
if param.includenuclear || param.includecell
    if param.includenuclear
       
        for i=1:length(nuclearparams)
            for j=1:n
                y2(j,i)=eval(['models{j}.nuclearShapeModel.' nuclearparams{i}]');
            end
        end

    end

    tparam = [];
    tparam.generatemeanshape = true;
    tparam.display = false;
    tparam.debug = false;

    havealready = length(nuclearparams);
    for j=1:n
        model = models{j};
        tparam.resolution.cell = model.cellShapeModel.resolution;
        tparam.resolution.objects = eval(['model.' protfield '.resolution']);
        % generate average nuclear instance
        
        if is3D
            instance = tp_genspsurf(model.nuclearShapeModel, tparam);
            tp_gennucshape_result = tp_gennucshape(instance, tparam);
            nucimg = tp_gennucshape_result.nucimg;
            nucsurf = tp_gennucshape_result.nucsurf;

            nucleus.nucimgsize = size(nucimg);
            nucleus.nucsurf = nucsurf; clear nucsurf;
            nucleus.nucimg = nucimg; clear nucimg;
            % generate cell shape instance from it
            ml_gencellshape3d_result = ml_gencellshape3d( ...
                model.cellShapeModel, nucleus, tparam );
            cellimg = ml_gencellshape3d_result.cellimg;
            nucimg = ml_gencellshape3d_result.nucimg;

            if param.includenuclear
                stats=rm_get3Dstats(nucimg);
                statnames = fieldnames(stats);
                for k=1:length(statnames)
                    y2(j, k+havealready) = eval(['stats.' statnames{k}]);
                    nuclearparams{k+havealready} = ['Nuclear ' statnames{k}];
                end

                clear nucimg;
            end

            if param.includecell
                stats=rm_get3Dstats(cellimg);
                statnames = fieldnames(stats);
                for k=1:length(statnames)
                    y3(j, k) = eval(['stats.' statnames{k}]);
                end
                
                cellparams = {};
                for k=1:length(statnames)
                   cellparams{end+1} = ['Cell ' statnames{k}];
                end
                %cellparams=statnames';

                clear cellimg;
            end
        
        end
    end
end

paramnames = [proteinparams nuclearparams cellparams];
yall = [y y2 y3];

%
% plot all of the parameters for all of the models individually, and plot
% the ones that change the most together
%

%save testing

%end

%%%%%%%%%%%%%%%%%%%%%% MODEL IDENTIFIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, '%%%%%% Model Identifiers\n');
for i=1:n
    identifiers = {'filename' 'name' 'id'};
    for j=1:3
        identifier = identifiers{j};
        model = models{i};
        if isfield(model, identifier)
            identifier_value = eval(sprintf('model.%s', identifier));
            fprintf( fileID, 'model%d.%s = ''%s'';\n', i, identifier, identifier_value );
        else
            fprintf( fileID, 'model%d.%s = ''%s'';\n', i, identifier, 'UNSET');
        end
    end
    fprintf( fileID, '\n' );
end

%%%%%%%%%%%%%%%%%%%%%% FIGURES 1 & 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model1_name = models{1}.name;
model2_name = models{2}.name;
fig_count = 0;
fig_names = {sprintf('Protein Model Parameters'), ...
    sprintf('Nucleus and Cell Model Parameters'), ...
    'Coefficient of Variation', ...
    'Fractional Distance', ...
    'Number of Objects'};
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
    I = imread( sprintf('image%d.png', fig_count) );
    I = imresize( I, 0.5 );
    imwrite( I, sprintf('image%d_thumbnail.png', fig_count) );
    fprintf( fileID, '%%%%%% %s\n', fig_names{fig_count});
    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="image%d.png" ><img src="image%d_thumbnail.png" /></a>\n', fig_count, fig_count );
    fprintf( fileID, '%% <figcaption><strong>Figure %d.</strong> model1 (vesicle1), model2 (vesicle2)</figcaption>\n', fig_count );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
end

%%%%%%%%%%%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures = [figures, figure];
fig_count = fig_count + 1;
normy = yall./repmat(mean(yall),n,1);
cv=std(normy);
[sorted,order]=sort(cv,'descend');
% find the ones whos CV is greater than 50%
idx = find(sorted>0.5);
mintodisplay = 3;
if length(idx) > mintodisplay
    order = order(sorted>0.5);
%    order = order(end:-1:1);
else
    order = order(1:mintodisplay);
%    order = order(end:-1:end-mintodisplay+1);
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
I = imread( sprintf('image%d.png', fig_count) );
I = imresize( I, 0.5 );
imwrite( I, sprintf('image%d_thumbnail.png', fig_count) );
fprintf( fileID, '%%%%%% %s\n', fig_names{fig_count});
fprintf( fileID, '%%\n' );
fprintf( fileID, '%% <html>\n' );
fprintf( fileID, '%% <figure>\n' );
fprintf( fileID, '%% <a href="image%d.png" ><img src="image%d_thumbnail.png" /></a>\n', fig_count, fig_count );
fprintf( fileID, '%% <figcaption><strong>Figure %d.</strong> model1 (vesicle1), model2 (vesicle2)</figcaption>\n', fig_count );
fprintf( fileID, '%% </figure>\n' );
fprintf( fileID, '%% </html>\n' );
fprintf( fileID, '%%\n' );

%%%%%%%%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.includeprot
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
I = imread( sprintf('image%d.png', fig_count) );
I = imresize( I, 0.5 );
imwrite( I, sprintf('image%d_thumbnail.png', fig_count) );
fprintf( fileID, '%%%%%% %s\n', fig_names{fig_count});
fprintf( fileID, '%%\n' );
fprintf( fileID, '%% <html>\n' );
fprintf( fileID, '%% <figure>\n' );
fprintf( fileID, '%% <a href="image%d.png" ><img src="image%d_thumbnail.png" /></a>\n', fig_count, fig_count );
fprintf( fileID, '%% <figcaption><strong>Figure %d.</strong> model1 (vesicle1), model2 (vesicle2)</figcaption>\n', fig_count );
fprintf( fileID, '%% </figure>\n' );
fprintf( fileID, '%% </html>\n' );
fprintf( fileID, '%%\n' );

%%%%%%%%%%%%%%%%%%%%% FIGURE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    I = imread( sprintf('image%d.png', fig_count) );
    I = imresize( I, 0.5 );
    imwrite( I, sprintf('image%d_thumbnail.png', fig_count) );
    fprintf( fileID, '%%%%%% %s\n', fig_names{fig_count});
    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="image%d.png" ><img src="image%d_thumbnail.png" /></a>\n', fig_count, fig_count );
    fprintf( fileID, '%% <figcaption><strong>Figure %d.</strong> model1 (vesicle1), model2 (vesicle2)</figcaption>\n', fig_count );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
end
end