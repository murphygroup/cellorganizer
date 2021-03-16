function feat = ml_gaborfeat(img,tmpdir)

%ML_GABORFEAT calculates gaobr features
%   FEAT = ML_GARBORFEAT(IMG) calculates 60 Gabor texture features 
%   for a 2D image IMG, which must be uint8.
%   FEAT = ML_GARBORFEAT(IMG,TMPDIR) allows user to specify a directory
%   TMPDIR for saving temporary files, which is '/tmp' for default.
 
% Copyright (C) 2006  Murphy Lab
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

%June, 2005 T. Zhao

if ~isa(img,'uint8')
    error('The input image must be uint8');
end

[tmppdir,tmpfile1]=fileparts(tempname);
[tmppdir,tmpfile2]=fileparts(tempname);

if ~exist('tmpdir','var')
    tmpdir=tempdir;
end

if( tmpdir(end) ~= filesep) 
    tmpdir(end+1) = filesep; 
end;

imgfile=[tmpdir tmpfile1];
featfile=[tmpdir tmpfile2];

fid=fopen(imgfile,'w');
fwrite(fid,img,'uchar');
fclose(fid);

bincmd=which('ml_gaborfeat');
slashpos=find(bincmd==filesep);
bincmd=[bincmd(1:slashpos(end-1)) 'bin' bincmd(slashpos(end):end-2)];

system([bincmd ' ' imgfile ' ' num2str(size(img,1)) ' '...
        num2str(size(img,2)) ' ' featfile]);
    
if ~exist(featfile,'file')
    feat = zeros(1,60)+NaN;
    return
end

fid=fopen(featfile,'r');
feat=fread(fid,60,'float');
fclose(fid);
feat=feat';

delete(imgfile);
delete(featfile);