function feat = tz_garborfeat(img,tmpdir)
%TZ_GARBORFEAT Obsolete.
%
%See also ML_GARBORFEAT

%function feat = tz_garborfeat(img,tmpdir)
%OVERVIEW
%   
%PARAMETERS
%   img - input image, must be uint8
%   tmpdir - 
%RETURN
%   feat - 
%DESCRIPTION
%   
%HISTORY
%   28-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('feat','tz_garborfeat','ml_garborfeat'));

if ~exist('tmpdir','var')
    tmpdir='/tmp';
end

imgfile=[tmpdir '/' 'img.gabortexture'];
featfile=[tmpdir '/' 'feat.gabortexture'];
fid=fopen(imgfile,'w');
fwrite(fid,img,'uchar');
fclose(fid);
gaborfeatcmd=which('tz_garborfeat');
unix([gaborfeatcmd(1:end-2) ' ' imgfile ' ' num2str(size(img,1)) ' '...
        num2str(size(img,2)) ' ' featfile]);
fid=fopen(featfile,'r');
feat=fread(fid,60,'float');
fclose(fid);
feat=feat';