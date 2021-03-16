function answer = check_if_files_exist_on_disk_and_link_them_LAT( string )
% CHECK_IF_FILES_EXIST_ON_DISK Simple helper function that takes a string
% from Galaxy and checks if these strings match local filenames

% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
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

delimiter = ',';
C = strsplit(string,delimiter);

if isempty(C)
    answer = false;
    return
end

answer = true;
for i=1:length(C)
    filename = C{i};
    if exist(filename)
        tiffname=[num2str(i) '.ome.tif'];
        disp(['File ' filename ' exists. Linking file to ' tiffname]);
        system(['ln -s ' filename ' ' tiffname]);
        % [status,ll]=system(['/usr/local/Cellar/bftools/showinf -omexml-only -nopix ' tiffname ' |grep "<Channel ID="']);
        [status,ll]=system(['showinf -omexml-only -nopix ' tiffname ' |grep "<Channel ID="']);
        % cmdout=strsplit(cmdout,'\n');
        % for ii=1:length(cmdout)
        % ll=cmdout;
        ll=strsplit(ll,'"');
        tif_path=ll{4};
        disp(tif_path);
        tt=strsplit(tif_path,'/');
        ttt=tt(1:end-1);
        tif_folder=strjoin(ttt,'/');
        tif_folder=strrep(tif_folder,' ','\ ');

        system(['mkdir -p ' tif_folder]);
        img=OME_loadchannel(tiffname,1);
        img2tif(img,tif_path);
        % end
    else
        warning(['File ' filename ' does not exist.']);
        answer = false;
    end
end

end
