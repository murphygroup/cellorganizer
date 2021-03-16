function tz_yeastfeats_p(datafile)
%TZ_YEASTFEATS_P Calculate features for yeast image.
%   TZ_YEASTFEATS_P(DATAFILE,RESULTDIR)
%   
%   See also

%   22-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 arguments are required')
end

addpath /home/tingz/matlab
tz_initpath

[fileid,filename] = ml_fopentemp('/home/tingz/tmp');

load(datafile)
nprocfiles = length(procfiles{1})
isbreak = 0;

prvfile = [prvfile filename];

while 1
    while unix(['mv ' ctrfile ' ' prvfile]) == 1
        pause(1)
    end
    
    while ~exist(prvfile,'file')
        pause(.1)
    end
    
    [fid,message] = fopen(prvfile,'r');
    n = fread(fid,1,'uint16');
    fclose(fid);
    
    if isempty(n) | length(n)==0
        isbreak = 1;
    end
    if n>nprocfiles
        isbreak =1;
    end
    
    if isbreak
        unix(['mv '  prvfile ' ' ctrfile]);
        break;
    end
    
    %%%%%%%%TO DO%%%%%%%%%%
    resultfile = [resultdir filesep 'feats' num2str(n) '.mat'];
    if ~exist(resultfile,'file');
%         features = n;
        img = ml_readimage(procfiles{3}{n});
        dnaimg = ml_readimage(procfiles{1}{n});
        [names, features, slfnames] = ml_featset( double(img), [], ...
            double(dnaimg), 'mcell_all',0,0,'quantile',...
            'rc');
        
        save(resultfile,'features');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    fid = fopen(prvfile,'w');
    n = n+1
    fwrite(fid,n,'uint16');
    fclose(fid)
    
    unix(['mv '  prvfile ' ' ctrfile]);
end