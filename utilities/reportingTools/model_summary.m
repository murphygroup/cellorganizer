function model_summary(model, varargin)
% model_summary prints useful information about a trained vesicle model to the
% console and saves plots of trained parameter distributions.
% 
% 
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% model                       trained vesicle model
% savedir                     directory to save plots
% filetype                    filetype for saved plots: type <help saveas> in the console for a list of options
% prefix                      prefix for savefiles, default is report
%
% Author: Tim Majarian (tmajarian@cmu.edu)
%
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


%% Setup
is3D = strcmpi(model.dimensionality,'3d');
if is3D
    compartment = model.proteinModel.cytonuclearflag;
else
    compartment = 'all';
end
if numel(varargin)==0
    filetype = 'pdf';
    savedir = pwd;
    prefix = 'report';
elseif numel(varargin)==1
    filetype = 'pdf';
    savedir = varargin{1};
    prefix = 'report';
elseif numel(varargin)==2
    filetype = varargin{2};
    savedir = varargin{1};
    prefix = 'report';
else 
    filetype = varargin{2};
    savedir = varargin{1};
    prefix = varargin{3};
end

if ~exist(savedir,'dir')
    try 
        mkdir(savedir)
    catch
        disp('Could not create save directory, aborting')
        return
    end
end


%% Print model info
try
    disp(sprintf('\n'))
    disp('Model Overview')
    disp(['Model Name: ' model.proteinModel.name])
    disp(['Model ID: ' model.proteinModel.id])
    disp(['Model filename: ' model.filename])
    disp(['Cell Line: ' model.cellline])
    disp(['Training info: ' model.description])
    disp(['Trained in CellOrganizer version ' model.version])
    disp(sprintf('\n'))
    disp('Model Specifications')
    disp(['Model type: ' model.proteinModel.type])
    disp(['Model class: ' model.proteinModel.class])
    disp(['Dimensionality: ' model.proteinModel.dimensionality]) 
    disp(['Training resolution: ' model.proteinModel.resolution])
    disp(['Training downsampling: ' model.info.downsampling_vector])
    disp([model.documentation.description ' on ' model.documentation.date])
    disp(['Author: ' model.documentation.author])
    disp(['Questions? Concerns? Direct to: ' model.documentation.mailing_list])
catch
    disp(sprintf('\n'))
    disp('Model Overview')
    disp(['Model Name: ' model.proteinModel.name])
    %disp(['Model ID: ' model.proteinModel.id])
    %disp(['Model filename: ' model.filename])
    disp(['Cell Line: ' model.cellline])
    %disp(['Training info: ' model.description])
    disp(['Trained in CellOrganizer version ' model.version])
    disp(sprintf('\n'))
    disp('Model Specifications')
    disp(['Model type: ' model.proteinModel.type])
    disp(['Model class: ' model.proteinModel.class])
    %disp(['Dimensionality: ' model.proteinModel.dimensionality]) 
    disp(['Training resolution: ' model.proteinModel.resolution])
    %disp(['Training downsampling: ' model.info.downsampling_vector])
    %disp([model.documentation.description ' on ' model.documentation.date])
    disp(['Author: ' model.documentation.author])
    %disp(['Questions? Concerns? Direct to: ' model.documentation.mailing_list])
end


%% Plot 2D param
if ~is3D
    %% plot frq distro
    objalpha = model.proteinModel.objectModel.numStatModel.alpha;
    objbeta = model.proteinModel.objectModel.numStatModel.beta;
    xrange = 0:1000;
    gampd = gampdf(xrange, objalpha, objbeta);

    figure;
    plot(xrange, gampd)
    title('Distribution of object quantity')
    xlabel('Number of objects')
    ylabel('Relative probability density')        
    set(gca,...
    'YTickLabel','')
    saveas(gcf,[savedir filesep prefix '_frq_distro'],filetype);
    
    %% plot size distro
    mu = model.proteinModel.objectModel.intensStatModel.mu;
    sig = model.proteinModel.objectModel.intensStatModel.sigma;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in x')
    saveas(gcf,[savedir filesep prefix '_size_distro'],filetype);
    
else
%% Plot 3D param
    %% plot x size distro
    mu = model.proteinModel.size.x.mu;
    sig = model.proteinModel.size.x.sigma;
    x_vals = (0:.01:mu*10);
    y_vals = lognpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in x')
    saveas(gcf,[savedir filesep prefix '_x_size_distro'],filetype);


    %% plot y size distro 1
    mu = model.proteinModel.size.y_x.a1;
    sig = model.proteinModel.size.y_x.b1;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>.0001);
    y_vals = y_vals(y_vals>.0001);
    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in y1')
    saveas(gcf,[savedir filesep prefix '_y1_size_distro'],filetype);


    %% plot y size distro 2
    mu = model.proteinModel.size.y_x.a2;
    sig = model.proteinModel.size.y_x.b2;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in y2')
    saveas(gcf,[savedir filesep prefix '_y2_size_distro'],filetype);


    %% plot z size distro 1
    mu = model.proteinModel.size.z_x.a1;
    sig = model.proteinModel.size.z_x.b1;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in z1')
    saveas(gcf,[savedir filesep prefix '_z1_size_distro'],filetype);


    %% plot z size distro 2
    mu = model.proteinModel.size.z_x.a2;
    sig = model.proteinModel.size.z_x.b2;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object size in z2')
    saveas(gcf,[savedir filesep prefix '_z2_size_distro'],filetype);


    %% Frequency plot
    mu = model.proteinModel.frequency.mu;
    sig = model.proteinModel.frequency.sigma;
    x_vals = (0:.01:mu*10);
    y_vals = normpdf(x_vals,mu,sig);
    x_vals = x_vals(y_vals>0.0001);
    y_vals = y_vals(y_vals>0.0001);

    figure;
    plot(x_vals,y_vals)
    title('Distribution of object frequency')
    saveas(gcf,[savedir filesep prefix '_obj_frq_distro'],filetype);

end
%% Position plots
if is3D
    [positionBeta, fractdist] = beta2posMap(model.proteinModel.position.beta);
else
    [positionBeta, fractdist] = beta2posMap(model.proteinModel.positionModel.beta);
end
        
objpdf = squeeze(sum(sum(positionBeta,1),3));
    
if strcmpi(compartment, 'cyto')
    domain = fractdist >= 0;
    objpdf(~domain) = [];
    objpdf = [zeros(1,1+sum(~domain)), objpdf ] ;
    fractdists = [fractdist(~domain) 0 fractdist(domain)];

elseif strcmpi(compartment, 'nuc')
    domain = fractdist <= 0;
    objpdf(~domain) = [];
    objpdf = [objpdf zeros(1,1+sum(~domain))];
    fractdists = [fractdist(domain) 0 fractdist(~domain)];

else
    fractdists = fractdist;
end
        
objpdf = objpdf./sum(objpdf);

figure;
plot(fractdists, objpdf, 'LineWidth',2)
title('Fractional distance between nucleus and membrane')
xlim([-.1,fractdists(size(objpdf,2))+.1]')
saveas(gcf,[savedir filesep prefix '_fract_dist'],filetype);
