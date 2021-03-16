function  script_kde_trans_model_alon(model, modelname, object_info, parentdir, param_in)

% Copyright (C) 2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

if ~exist('param','var')
    param_in = [];
end

param_in = ml_initparam(param_in, struct('skipsynth', false, ...
    'keepinds', true(size(model.cellShapeModel.positions,1),1), ...
    'suffix', ''));


figuredir = [parentdir filesep 'figures'];
if ~exist(figuredir, 'dir')
    mkdir(figuredir)
end

tempdir = [figuredir filesep mfilename '_' modelname(1:end-4) param_in.suffix];
if ~exist(tempdir, 'dir')
    mkdir(tempdir)
end

syntrackdir = [tempdir filesep 'syntracks'];
if ~exist(syntrackdir, 'dir')
    mkdir(syntrackdir)
end

framedir = [tempdir filesep 'frames'];
if ~exist(framedir, 'dir')
    mkdir(framedir)
end

kerneldir = [tempdir filesep 'kernel'];
if ~exist(kerneldir, 'dir')
    mkdir(kerneldir)
end

pos_no_nan = model.cellShapeModel.positions;

pos_orig = model.cellShapeModel.positions;
pos = model.cellShapeModel.positions;

pos(~param_in.keepinds,:) = nan;
pos_orig(~param_in.keepinds,:) = nan;

uconditions = unique(object_info.obj_condition);
ndat = size(pos,1);
%remove outliers

model.cellShapeModel.positions = pos;

shapespacefile = [tempdir filesep 'shapeSpace.tiff'];
shapespacefile_param = [shapespacefile '.mat'];
if ~exist(shapespacefile, 'file')
    
    param.subsize = 2400;
    figure,
    
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [3.66, 3.66])
    
    rng(0)
    param.synthorder = randperm(model.cellShapeModel.numimgs);
    
    param.rotate = 45;
    
    param = showShapeSpaceFigure(model, object_info.obj_condition, param);
    
    %         set(get(gcf, 'CurrentAxes'), 'FontSize', 8)
    %         set(get(gcf, 'CurrentAxes'),'position',[0.05 0.05 .95 .95],'units','normalized')
    %     %     set(get(h(i), 'CurrentAxes'), 'Axes', 'tight')
    %
    
    set(get(gcf, 'CurrentAxes'),'position',[0.08, 0.08, .9, 0.9],'units','normalized')
    set(get(gcf, 'CurrentAxes'), 'FontSize', 6)
    print(shapespacefile, gcf, '-dtiff', '-r600')
    
    %         saveas(gcf, shapespacefile, 'tif')
    close(gcf)
    
    save(shapespacefile_param, 'param', '-v7.3')
else
    load(shapespacefile_param)
end

[dpos_grid] = 10.^(-6:0.05:2);

pos_min = min(pos);
pos_max = max(pos);

[meshx, meshy] = meshgrid(linspace(-0.5, 0.5, 300), linspace(-0.5, 0.5, 300));

all_pos = [meshx(:), meshy(:), zeros(numel(meshx), size(pos,2)-2) ];

shapespace_contour_file = [ tempdir filesep 'shapeSpace_contour.tif'];
densities_file = [tempdir filesep 'densities.mat'];

densities = zeros(size(model.cellShapeModel.positions,1),1);
if ~exist(shapespace_contour_file, 'file')
    figure,
    set(get(gcf, 'CurrentAxes'), 'FontSize', 6)
    %         set(get(gcf, 'CurrentAxes'),'position',[0,0,1,1],'units','normalized')
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [3.66, 3.66])
    for condition_ind = 1:length(uconditions)
        ucondinds = (object_info.obj_condition == uconditions(condition_ind)) & param_in.keepinds' ;
        pos_cond = model.cellShapeModel.positions(ucondinds,:);
        
        dmat_cond = squareform(pdist(pos_cond));
        
        dmat_cond(logical(eye(size(dmat_cond)))) = inf;
        
        k_best = find_kernel(dmat_cond, dpos_grid);
        
        dists = squareform(pdist(pos_cond));
        w = zeros(size(dists));
        w(:) = mvnpdf(dists(:), 0, k_best);
        w = sum(w,2);
        
        densities(ucondinds) = w;
        
        bound_50 = prctile(w, 50);
        bound_90 = prctile(w, 90);
        
        all_dists = pdist2(all_pos, pos_cond);
        
        w_dists = zeros(size(all_dists));
        w_dists(:) = mvnpdf(all_dists(:), 0, k_best);
        w_dists = sum(w_dists,2);
        
        w_dists = reshape(w_dists, size(meshx));
        
        %              scatter(pos_cond(:,1), pos_cond(:,2), 50, w, '.')
        %              hold on
        contour(meshx, meshy, w_dists, bound_50, 'LineColor', param.colors(condition_ind,:))
        hold on
        contour(meshx, meshy, w_dists, bound_90,'LineWidth',1.5, 'LineColor', param.colors(condition_ind,:));
    end
    axis equal
    axis(param.axis)
    set(get(gcf, 'CurrentAxes'),'position',[0.08, 0.08, .9, 0.9],'units','normalized')
    set(get(gcf, 'CurrentAxes'), 'FontSize', 6)
    
    save(densities_file, 'densities')
    %         saveas(gcf, shapespace_contour_file, 'tif')
    print(shapespace_contour_file, gcf, '-dtiff', '-r600')
    close(gcf)
    
end


naninds = find(any(isnan(pos),2));
%for each tagged protein
for condition_ind = 1:length(uconditions)
    
    
    %find the objects that correspond to this condition
    ucondinds = object_info.obj_condition == uconditions(condition_ind) ;
    objinds = find(ismember(object_info.transition_obj_inds(:,1), find(ucondinds)));
    
    objinds = objinds(all(~ismember(object_info.transition_obj_inds(objinds,:), naninds),2));
    
    pos0 = pos(object_info.transition_obj_inds(objinds,1),:);
    pos1 = pos(object_info.transition_obj_inds(objinds,2),:);
    
    inds_all{condition_ind} = object_info.transition_obj_inds(objinds,1);
    
    movieids = object_info.umovieid(object_info.transition_obj_inds(objinds,1));
    
    
    dmat_embed = squareform(pdist(pos0));
    
    kernel_sub_dir = [kerneldir filesep 'k' num2str(condition_ind)];
    if ~exist(kernel_sub_dir, 'dir')
        mkdir(kernel_sub_dir)
    end
    
    
    kernelfile = [kernel_sub_dir filesep 'kernel.mat'];
    if ~exist(kernelfile, 'file')

        fvals = nan(1, length(dpos_grid));
        
        for i = 1:length(dpos_grid(:))
            
            kernelfile_search = [kernel_sub_dir filesep 'k_' num2str(i)];
            
            if ~exist([kernelfile_search, '.mat'], 'file')
                [can_start, ~, final_exists, temp_name] = chunk_start( kernelfile_search, '.mat');
                if can_start && ~final_exists
                    fval = opt_func(pos0, pos1, dmat_embed, [dpos_grid(i)], true);
                    fvals(i) = fval;
                    
                    save([kernelfile_search '.mat'], 'fval')
                end
                
            else
                load([kernelfile_search '.mat'])
                fvals(i) = fval;
            end
        end
        
        if any(isnan(fvals))
            return;
        end
        
        
        [~, ind] = min(fvals);
        
        
        kernelfile
        
        disp('Computing kernel');
        [ksize, fval, exitflag] = fminsearch(@(x) opt_func(pos0, pos1, dmat_embed, x, movieids), dpos_grid(ind), optimset('Display','iter', 'MaxIter', 2000));
        save(kernelfile, 'ksize', 'fval', 'exitflag', 'fvals', 'dpos_grid')
        
    else
        load(kernelfile)
    end
    
    kernels(condition_ind) = ksize;
    fvals_all{condition_ind} = fvals;
    
    pos0_all{condition_ind} = pos0;
    pos1_all{condition_ind} = pos1;
end

kernel_response_file = [tempdir filesep 'kernel_err.png'];
if ~exist(kernel_response_file, 'file')
    colors = lines(length(fvals_all));
    figure('color', 'w'),
    hold on
    for i = 1:length(fvals_all)
        plot(log10(dpos_grid), fvals_all{i}, 'color', colors(i,:));
    end
    
    legend(gca, cellstr(num2str([1:length(fvals_all)]')))
    
    xlabel('log10(kernel size)')
    ylabel('HOO MSE')
    
    saveas(gcf, kernel_response_file, 'png')
end
% ksize = [ksize(1), 0;0, ksize(2)];


%Print off the vector field images
dimrange = [1:2];
pts_per_dim = 30;

p_low = prctile(pos, 2);
p_high = prctile(pos, 98);

interval = min((p_high(dimrange) - p_low(dimrange))/pts_per_dim);

numz = 5;
zbound = max(abs([prctile(pos(:,3), 10), prctile(pos(:,3), 90)]));
zspace = linspace(-zbound, zbound, numz);


[xpos,ypos] = meshgrid([p_low(1):interval:p_high(1)], ...
    [p_low(2):interval:p_high(2)]); ...
    

for z = 1:length(zspace)
    coords = [xpos(:), ypos(:), ones(size(ypos(:))).*zspace(z), zeros(length(xpos(:)), size(pos,2)-3)];
    
    %get the predicted displacements over a range of coordinates
    for condition_ind = 1:length(uconditions)
        
        pos0 = pos0_all{condition_ind};
        pos1 = pos1_all{condition_ind};
        
        coords1 = zeros(size(coords));
        
        for i = 1:length(coords)
            w = mvnpdf([pdist2(coords(i,:), pos0)'], 0, kernels(condition_ind));
            w = w./sum(w(:));
            coords1(i,:) = pos1'*w;
        end
        
        delta{condition_ind}{z}= coords1 - coords;
        velocity{condition_ind}{z} = sqrt(sum((coords1 - coords).^2,2));
    end
end

deltas = [delta{:}];

deltas = vertcat(deltas{:});
d_low = prctile(deltas(:), 2);
d_high = prctile(deltas(:), 98);

d_high = max(abs([d_low, d_high]));
d_low = -d_high;


velocities = [velocity{:}];

velocities = vertcat(velocities{:});
v_low = prctile(velocities(:), 2);
v_high = prctile(velocities(:), 98);

%     v_high = max(abs([v_low, v_high]));
v_low = 0;

vecfielddir = [tempdir filesep 'vecfields'];
if ~exist(vecfielddir, 'dir')
    mkdir(vecfielddir)
end

%for each condition, print off a vector field image
for condition_ind = 1:length(uconditions)
    for z = 1:length(zspace)
        vecfieldimg = [vecfielddir filesep 'vecfield_cond' num2str(condition_ind) '_z' num2str(z) '.tif'];
        vecfieldimg_velocity = [vecfielddir filesep 'vecfield_velocity_cond' num2str(condition_ind) '_z' num2str(z) '.tif'];
        if ~exist(vecfieldimg, 'file')
            
            pos0 = pos0_all{condition_ind};
            pos1 = pos1_all{condition_ind};
            
            d = zeros(size(xpos));
            d(:) = delta{condition_ind}{z}(:,3);
            
            figure('color', 'w', 'PaperUnits', 'Inches', 'PaperSize', [2.30, 4] )
            [~,h] = contourf(xpos, ypos, d, linspace(d_low, d_high, 1000));
            
            set(h,'LineStyle','none');
            
            colormap(othercolor('BuDRd_18',1000))
            
            hold on
            h = quiver(coords(:,1), coords(:,2), delta{condition_ind}{z}(:,1), delta{condition_ind}{z}(:,2), 3, 'k');
            
            
            set(h, 'LineWidth', 0.25)
            
            %                 xlabel('MDS 1')
            %                 ylabel('MDS 2')
            grid off
            axis equal
            axis tight
            
            hold on
            scatter(pos0(:,1), pos0(:,2), 50, 'k', '.')
            
            title(['Condition ' num2str(condition_ind) ' z = ' num2str(zspace(z))]);
            
            set(gca, 'FontSize', 6)
            set(get(gcf, 'CurrentAxes'),'position',[0.08, 0.08, 0.9, 0.9],'units','normalized')
            set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [2.30, 5] )
            
            print(vecfieldimg, gcf, '-dtiff', '-r600')
            %                 saveas(gcf, vecfieldimg, 'tif')
            close gcf;
            
            v = zeros(size(xpos));
            v(:) = velocity{condition_ind}{z};
            
            figure('color', 'w', 'PaperUnits', 'Inches', 'PaperSize', [2.30, 4] )
            [~,h] = contourf(xpos, ypos, v, linspace(v_low, v_high, 1000));
            
            set(h,'LineStyle','none');
            
            colormap(othercolor('OrRd9',1000))
            
            hold on
            h = quiver(coords(:,1), coords(:,2), delta{condition_ind}{z}(:,1), delta{condition_ind}{z}(:,2), 3, 'k');
            
            
            set(h, 'LineWidth', 0.25)
            
            %                 xlabel('MDS 1')
            %                 ylabel('MDS 2')
            grid off
            axis equal
            axis tight
            
            hold on
            scatter(pos0(:,1), pos0(:,2), 50, 'k', '.')
            
            title(['Condition ' num2str(condition_ind) ' z = ' num2str(zspace(z))]);
            
            set(gca, 'FontSize', 6)
            set(get(gcf, 'CurrentAxes'),'position',[0.08, 0.08, 0.9, 0.9],'units','normalized')
            set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [2.30, 5] )
            
            print(vecfieldimg_velocity, gcf, '-dtiff', '-r600')
            %                 saveas(gcf, vecfieldimg, 'tif')
            close gcf;
            
        end
    end
end

h = figure('PaperUnits', 'Inches', 'PaperSize', [2.30,1]);
caxis([d_low, d_high]);
colormap(othercolor('BuDRd_18',1000))
colorbar
hold on
set(gca, 'FontSize', 6)
plot(0,0)
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [2.30,1])

print([vecfielddir filesep 'colorbar.tif'], h, '-dtiff', '-r600')


h = figure('PaperUnits', 'Inches', 'PaperSize', [2.30,1]);
caxis([v_low, v_high]);
colormap(othercolor('OrRd9',1000))
colorbar
hold on
set(gca, 'FontSize', 6)
plot(0,0)
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [2.30,1])

print([vecfielddir filesep 'colorbar_velocity.tif'], h, '-dtiff', '-r600')



nframes = 69;
numsyntracks = 5;


model.nuclearShapeModel = model.cellShapeModel;

convhulldatfile = [tempdir filesep 'convhulldat.mat'];

if ~exist(convhulldatfile, 'file')
    param_embed.positions = pos_orig;
    [ embedding_model ] = embed_distance_matrix( pos_orig, param_embed);
    
    save(convhulldatfile, 'embedding_model')
else
    load(convhulldatfile)
end

convhull.inds = embedding_model.convex_hull;
convhull.tes = embedding_model.tessellation;
convhull.positions = embedding_model.positions;


model.cellShapeModel.positions = pos_no_nan;
model.cellShapeModel.tessellation = convhull.tes;
model.cellShapeModel.convex_hull = convhull.inds;

model.nuclearShapeModel = model.cellShapeModel;



if ~param_in.skipsynth
    for condition_ind = 1:length(uconditions)
        for i = 1:numsyntracks
            ind = sub2ind([length(uconditions), numsyntracks], condition_ind, i);
            rng(ind);
            
            trackfile = [syntrackdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i)];
            
            if ~exist([trackfile '.mat'], 'file')
                [can_start, ~, final_exists, temp_name] = chunk_start( trackfile, '.mat');
                
                if can_start & ~final_exists
                    disp(['Generating tracks: ' trackfile]);
                    
                    pos0 = pos0_all{condition_ind};
                    pos1 = pos1_all{condition_ind};
                    
                    mpos = mean([pos0;pos1]);
                    cpos = cov([pos0;pos1]);
                    
                    fpos = nan(nframes,size(pos0,2));
                    fpos(1,:) = mvnrnd(mpos, cpos, 1);
                    a = nan;
                    while isnan(a)
                        fpos(1,:) = mvnrnd(mpos, cpos, 1);
                        a = tsearchn(convhull.positions,convhull.tes,fpos(1,:));
                    end
                    
                    for j = 2:nframes
                        disp([num2str(j) filesep num2str(nframes)]);
                        %compute the weights
                        d_dist = pdist2(fpos(j-1,:), pos0)';
                        
                        w = mvnpdf(d_dist,0, ksize);
                        w = w./sum(w);
                        
                        %average displacement
                        m_disp = (pos1'*w)';
                        
                        residuals = pos1 - repmat(m_disp, [size(pos1,1),1]);
                        % Weighted Covariance Matrix
                        cov_disp = residuals' * (residuals .* repmat(w, 1, size(residuals,2)));
                        cov_disp = 0.5 * (cov_disp + cov_disp');
                        
                        fpos(j,:) = m_disp;
                        
                        a = nan;
                        while isnan(a)
                            fpos(j,:) = mvnrnd(m_disp, cov_disp);
                            a = tsearchn(convhull.positions,convhull.tes,fpos(j,:));
                        end
                        
                    end
                    save(trackfile, 'fpos')
                    chunk_finish(temp_name)
                end
            end
        end
    end
end

trackplotdir = [tempdir filesep 'trackplots'];

if ~exist(trackplotdir, 'dir')
    mkdir(trackplotdir)
end

colors = jet(length(uconditions));

randinds = randperm(size(pos,1));
for condition_ind = 1:length(uconditions)
    for i = 1:numsyntracks
        ind = sub2ind([length(uconditions), numsyntracks], condition_ind, i);
        rng(ind);
        
        trackfile = [syntrackdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i) '.mat'];
        trackplotfile = [trackplotdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i) '.tif'];
        trackplot_save_dir = [trackplotdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i)];
        if ~exist(trackplot_save_dir, 'dir')
            mkdir(trackplot_save_dir)
        end
        
        if exist(trackfile, 'file') && ~exist(trackplotfile, 'file') && ~param_in.skipsynth
            load(trackfile)
            figure('color', 'w')
            scatter(pos(:,1), pos(:,2), '.')
            hold on,
            plot(fpos(:,1), fpos(:,2), 'k')
            
            saveas(gcf, trackplotfile, 'tif')
            title(['Condition ' num2str(condition_ind) ', track ' num2str(i)])
            xlabel('mds1')
            ylabel('mds2')
            close gcf
            
            
            figure,
            scatter(pos(randinds,1), pos(randinds,2), 50, colors(object_info.obj_condition(randinds),:), '.' )
            axis([-0.3, 0.3, -0.15,0.15])
            hold on
            for j = 1:size(fpos,1)
                if j > 1
                    plot(fpos((j-1):j,1), fpos((j-1):j,2),'k', 'LineWidth',1)
                end
                scatter(fpos(j,1), fpos(j,2), 100, '.k')
                
                saveas(gcf, [trackplot_save_dir filesep 'frame' num2str(j) '.tif'], 'tif')
            end
        end
        
    end
end

interp_intervals = [0.5, 0.2];

for interval = 1:length(interp_intervals)
    framedir = [tempdir filesep 'frames_dir' num2str(interval) filesep];
    if ~exist(framedir, 'dir')
        mkdir(framedir)
    end
    
    moviedir = [tempdir filesep 'movies_' num2str(interval) filesep];
    
    if ~exist(moviedir, 'dir')
        mkdir(moviedir)
    end
    
    
    if ~param_in.skipsynth
        for condition_ind = 1:length(uconditions)
            for i = 1:numsyntracks
                trackfile = [syntrackdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i) '.mat'];
                
                if exist(trackfile, 'file')
                    load(trackfile)
                    
                    fpos = interp1(fpos, 1:interp_intervals(interval):size(fpos,1));
                    
                    for j = 1:size(fpos,1)
                        framefile = [framedir filesep 'cond' num2str(condition_ind) '_f' num2str(i) '_' num2str(j) '.mat'];
                        
                        tmpfile = [framefile '.tmp'];
                        
                        if ~exist(framefile, 'file') && ~exist(tmpfile, 'file')
                            system(['touch "' tmpfile '"'])
                            disp(framefile)
                            
                            modelparam.position = fpos(j,:);
                            modelparam.framefolder = [framedir filesep 'frame' num2str(condition_ind) '_f' num2str(i) '_' num2str(j)];
                            
                            system(['rm ' modelparam.framefolder filesep 'temp' filesep '*.tmp'])
                            
                            img = getfield(model2img({model}, modelparam), 'imgs');
                            
                            save(framefile, 'img', 'param')
                            
                            system(['rm ' tmpfile])
                            system(['rm -r ' modelparam.framefolder])
                            
                        end
                    end
                end
            end
        end
    end
    
    imsize = size(model.cellShapeModel.imfunc(1));
    missingimgs = cell(numsyntracks,1);
    
    frame_panels = [4, 6];
    
    for condition_ind = 1:length(uconditions)
        for i = 1:numsyntracks
            trackfile = [syntrackdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i) '.mat'];
            if exist(trackfile, 'file')
                load(trackfile)
                
                fpos = interp1(fpos, 1:interp_intervals(interval):size(fpos,1));
                
                missingimgs{i} = false(1, size(fpos,1));
                
                imgs = cell(1, size(fpos,1));
                
                trackplot_save_dir = [trackplotdir filesep 'syntracks' num2str(condition_ind) '_' num2str(i)];
                for j = 1:size(fpos,1)
                    
                    framefile = [framedir filesep 'cond' num2str(condition_ind) '_f' num2str(i) '_' num2str(j) '.mat'];
                    
                    if exist(framefile, 'file')
                        try
                            load(framefile)
                        catch
                            img = [];
                        end
                        
                        if ~isempty(img)
                            img1 = double(sum(img{1},3)>0).*1;
                            img2 = double(sum(img{2},3)>0).*0.5;
                            
                            imgs{j} = img1+img2;
                            %                         imgs{j} = ml_findmainobj(double((img{1}+img{2})>0));
                        else
                            system(['rm ' framefile])
                            imgs{j} = zeros(imsize);
                            missingimgs{i}(j) = true;
                        end
                    else
                        imgs{j} = zeros(imsize);
                        missingimgs{i}(j) = true;
                    end
                    
                    plotfile = [trackplot_save_dir filesep 'frame' num2str(j) '.tif'];
                    if exist(plotfile, 'file')
                        trackplotimg = ml_readimage(plotfile);
                        trackplotimgs{j} = [imresize(repmat(imgs{j},[1,1,3]).*255, size(imgs{j}).*2, 'nearest') imresize(trackplotimg, size(imgs{j}).*2)];
                    else
                        trackplotimgs{j} = [repmat(imresize(repmat(imgs{j},[1,1,3]).*255, size(imgs{j}).*2, 'nearest'), [1,2,1])];
                    end
                end
                
                try
                    [~, croprange] = cropimg(sum(cat(3,imgs{:}),3), 5);
                    imgs_fig = cellfun(@(x) x(croprange(1):croprange(2), croprange(3):croprange(4)), imgs, 'UniformOutput', false);
                    
                    imgs_fig_final = cell(frame_panels);
                    imgs_fig_final(:) = imgs_fig(1:prod(frame_panels));
                    imgs_fig_final = imgs_fig_final';
                    
                    imgs_rows = cell(size(imgs_fig_final,1),1);
                    for k = 1:size(imgs_fig_final,1)
                        imgs_rows{k} = [imgs_fig_final{k,:}];
                    end
                    
                    imgs_fig_final = vertcat(imgs_rows{:});
                    
                    
                    imwrite(imgs_fig_final.*1.25,  [moviedir filesep 'movie_frames' num2str(condition_ind) '_' num2str(i) '.tif'], 'tif')
                    
                    imgs2movie(imgs, [moviedir filesep 'movie_' num2str(condition_ind) '_' num2str(i) '.avi'], 10)
                catch
                    warning('Unable to save plot.');
                end
                
                try
                    imgs2movie(trackplotimgs, [moviedir filesep 'movie_plot_' num2str(condition_ind) '_' num2str(i) '.avi'], 10)
                catch
                    warning('Unable to save plot.');
                end
            end
        end
    end
end
end

function [diff, pred_pos1] = opt_func(pos0,pos1,dmat,h, cv)


if ~exist('dmat', 'var') || isempty(dmat)
    dmat = squareform(pdist(pos0));
end

if ~exist('cv', 'var') || isempty(cv)
    cv = true;
end

%if not positive definite, the kernel is no good
if h <= 0
    diff = inf;
    pred_pos1 = inf;
    return
end

pred_pos1 = inf(size(pos0));

for i = 1:size(pos0,1)
    %     t_dist(t_dist > 0.5) = t_dist(t_dist > 0.5)-1;
    %     t_dist(t_dist < -0.5) = t_dist(t_dist < -0.5)+1;
    %
    ddist = dmat(i,:)';
    
    w = mvnpdf(ddist, 0, h);
    
    if isnumeric(cv)
        w(cv == cv(i)) = 0;
    elseif cv
        w(i) = 0;
    end
    
    if all(w == 0)
        w(:) = 1;
    end
    
    w = w./sum(w(:));
    
    pred_pos1(i,:) = pos1'*w;
    
end

diff = mean(sum((pos1-pred_pos1).^2,2),1);
end

function k_best = find_kernel(dmat_embed, kernel_grid)
fvals = nan(1, length(kernel_grid));

for i = 1:length(kernel_grid(:))
    
    fval = opt_func2(dmat_embed, kernel_grid(i));
    fvals(i) = fval;
    
end

[~, ind] = max(fvals);

%     disp('Computing kernel');
[k_best, fval, exitflag] = fminsearch(@(x) -opt_func2(dmat_embed, x), ...
    kernel_grid(ind), optimset('Display','off', 'MaxIter', 200));
end

function [log_likelihood] = opt_func2(dists,h)

%if not positive definite, the kernel is no good
if h <= 0
    log_likelihood = 0;
    
    return
end

log_likelihood = kernel_pred(dists, h);
end

function log_likelihood = kernel_pred(dists,h)
w = zeros(size(dists));
w(:) = mvnpdf(dists(:), 0, h);

log_likelihood = sum(log10(sum(w,2)),1);
end