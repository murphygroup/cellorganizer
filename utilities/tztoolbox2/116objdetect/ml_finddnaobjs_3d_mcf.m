function ml_finddnaobjs_3d_mcf(dirname,savedir,loadthre,classes)
%ML_FINDDNAOBJS_3D_MCF Find nuclear objects in 3d image files from classes
%   ML_FINDIMGOBJS_3D_MCF(DIRNAME,SAVEDIR,LOADTHRE) finds object of the 3D 
%   images under the directiroy DIRNAME, which has a hierachical structure.
%   The first level subdirectories of DIRNAME are classes, followed by cells
%   and channels. LOADTHRE is the threshold do remove bright pixels. See
%   ML_FINDIMGOBJS_3D for more details.
%   
%   ML_FINDIMGOBJS_3D_MCF(DIRNAME,SAVEDIR,LOADTHRE,CLASSES) specify what
%   classes will be processed by the [string array] CLASSES, in which each
%   element is the name of a class.
%
%   See also

%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.  
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if ~exist('classes','var')
    %class_names=tz_cleandirs(mv_dir(dirname));
    class_names = tz_ls(dirname,'dir');
else
    class_names = classes;
end

if ~exist(savedir,'dir')
    unix(['mkdir ' savedir]);
end

script = ['ml_finddnaobjs_3d_mcf(''' dirname ''',''' savedir ''','];
if isempty(loadthre)
    script = [script '[])'];
else
    script = [script num2str(loadthre) ')'];
end

 
NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = class_names{class};
    classdir = [savedir filesep class_name];
    if ~exist(classdir,'dir')
        mkdir(classdir);
    end
    
    %status=mkdir(savedir,class_name);
    objectdir=[classdir '/dnaobj'];
    
    if ~exist('objectdir','dir')
        mkdir(objectdir);
    end
    
    %cells=tz_cleandirs(mv_dir([dirname '/' class_name]));
    cells = tz_ls([dirname '/' class_name '/cell*'],'dir');
    
    ncell=length(cells);
    for N=1:ncell
       % status=mkdir(saveclassdir,cells{N});
        %savecelldir=[saveclassdir '/' cells(N)];
        celldir=[dirname '/' class_name '/' cells{N}];
        protdir = [celldir '/dna'];
        cropfiles = tz_ls([celldir '/crop/*.tif*']);
        for k=1:length(cropfiles)
            savefile = [objectdir '/' cells{N} 'crop' num2str(k) '.mat'];
            if ~exist(savefile,'file')
                objects=ml_findimgobjs_3d(protdir,[celldir '/crop/' ...
                                    cropfiles{k}],loadthre);
                %save([saveclassdir '/' cells{N} '.mat'],'objects')
                tz_save(savefile,{objects},{'objects'},script,['3D ' ...
                                    'objects']);
            end
        end     
    end
end
