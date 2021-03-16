function cropconf = ...
    tz_majorslice_mcf(dirname,channel,savedir,ext,refchannel,isshow)
%TZ_MAJORSLICE_MCF Finding major slices in MCF 3D stacks.
%   CROPCONF = TZ_MAJORSLICE_MCF(DIRNAME,CHANNEL,SAVEDIR,EXT) extract a
%   slice with most intensities. See TZ_MAJORSLICE for details. This is a
%   batch processing function. All images with extension EXT under DIRNAME
%   (directly or indirectly) will be processed. The results will be saved 
%   in SAVEDIR. CHANNEL is the channel to be intensity calculation:  
%       'cell' - cell channel
%       'dna' - dna channel
%       'prot' - protein channel
%   
%   CROPCONF = TZ_MAJORSLICE_MCF(DIRNAME,CHANNEL,SAVEDIR,EXT,REFCHANNEL)
%   let user specify the reference channel for finding major slice. This
%   means that the major slice is searched in REFCHANNEL.The files under
%   CHANNEL and REFCHANNEL should have the same names.
%
%   CROPCONF = TZ_MAJORSLICE_MCF(DIRNAME,CHANNEL,SAVEDIR,EXT,REFCHANNEL,0)
%   turning displaying images off. 
%
%   See also

%   ??-???-???? Initial write T. Zhao
%   23-OCT-2004 Modified T. Zhao
%       - add channel parameter
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('4-6 arguments are required')
end

%find all class list
class_names=tz_cleandirs(ml_dir(dirname));

if ~exist(savedir,'dir')
    unix(['mkdir ' savedir]);
end

if ~exist('refchannel','var')
    refchannel = channel;
end

if ~exist('isshow','var')
    isshow = 1;
end

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;

cropconf={};
k=1;

for class = first_class : last_class
    class_name = class_names{class};
    if strmatch(class_name,{'rohrer','UCE'})
        continue;
    end
    status=mkdir(savedir,class_name);
    status=mkdir([savedir '/' class_name],channel);
    saveclassdir=[savedir '/' class_name '/' channel];
    
    cells=tz_cleandirs(ml_dir([dirname '/' class_name]),{'raw'});
    ncell=length(cells);
    for N=1:ncell
        celldir=[dirname '/' class_name '/' cells{N}];
        chandir = [celldir '/' channel];
        refchandir = [celldir '/' refchannel];
        savefile=[saveclassdir '/' cells{N} '_' channel '.mat'];
        cropfiles = ml_dir([celldir '/crop/*.tif*']);
        
        %record possible crop errors
        if length(cropfiles)>1
            cropconf{k}=celldir;
            k=k+1;
        end
        
        if ~exist(savefile,'file')
            %load mask
            mask=imread([celldir '/crop/' cropfiles{1}]);
            
            %find major slice
            [filename,selimg]=tz_findmajorslice(refchandir,ext,mask);
            if ~strcmp(chandir,refchandir)
                selimg = ml_readimage([chandir '/' filename]);
                selimg(mask==0)=0;
            end
            
            if(isshow==1)
                imshow(selimg,[]);
                title([chandir '/' filename])
                drawnow
            end
            
%             save(savefile,'selimg','filename');
            script = ['tz_majorslice_mcf(''' dirname ''',''' channel ...
                ''',''' savedir ''',''' ext ''',''' refchannel ...
                ''',' num2str(isshow) ')'];
            comments = 'major slice';
            tz_save(savefile,{selimg,filename},{'selimg','filename'}, ...
                script,comments);
        else
            tmp = load(savefile);
            imgfile = [chandir '/' tmp.filename];
            cropfile = [celldir '/crop/' cropfiles{1}];
            saveimgdir=[savedir '/' class_name '/' 'org' channel];
            if ~exist(saveimgdir,'dir')
                mkdir(saveimgdir);
            end
            copyfile(imgfile,[saveimgdir '/' cells{N} '.tif']);
            cropimgdir = [savedir '/' class_name '/crop'];
            if ~exist(cropimgdir,'dir')
                mkdir(cropimgdir);
            end
            copyfile(cropfile,[cropimgdir '/' cells{N} '.tif']);
        end
    end
end
