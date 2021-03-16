function answer = get_model_file_latest_version( model_file )
%GET_MODEL_FILE_LATEST_VERSION Retrieves latest version of a model file

% Copyright (C) 2007-2014 Murphy Lab
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

answer = false;

%get model files version list
cellorganizer_website = 'http://www.cellorganizer.org/Downloads/models';
model_files_version_list = [cellorganizer_website filesep 'models.txt'];
model_files_version_list_file = [ pwd filesep 'models.txt' ];
if exist( model_files_version_list_file );
    delete( model_files_version_list_file );
end

%get file containing the list of model files latest version
urlwrite( model_files_version_list, [ pwd filesep 'models.txt' ] );

%read the file contents
fid = fopen( model_files_version_list_file, 'rt');
lines = textscan(fid,'%[^\n]'); %reads line by line
lines = lines{1};
fclose(fid);

%generate superid from model file
superid = generate_superid_from_model_file( model_file );

%extract version from super id
version = superid(findstr( superid, 'version' )+8:end);
 
[ path, filename, extension ] = fileparts( model_file );

%build query string
query_string = superid(1:strfind(superid, 'version' )-2);

%search query string in list
if ~isempty( lines )
    for i=1:1:length(lines)
        line = lines{i};
        if findstr( query_string, line )
            %result from the query
            query_result = line;
            %extract version number
            query_result_version_number = query_result( findstr( query_result, 'version' )+8:end-4 );
            
            clear line
            
            break
        end
    end
end

%now that we have all the information we need to decide if we should
%upgrade this file or not
if upgrade( version, query_result_version_number )
    disp( 'Upgrading model file' );
    disp( 'Backing up older version of the file using their superid' );
    disp( [ 'Renaming ' model_file ' to ' superid ] );
    movefile( model_file, [ path filesep superid '.mat'] );
    disp( 'Retrieving and replacing with newer version' );
    urlwrite( [cellorganizer_website filesep query_result], model_file );
    answer=true;
else
    disp( 'Upgrade is not neccesary. Latest version found on disk.' );
end

if exist( model_files_version_list_file );
    delete( model_files_version_list_file );
end

end%get_model_file_latest_version

function answer = upgrade( version, latestVersion )
version = sscanf(version, '%d.%d.%d')';
if length(version) < 3
   version(3) = 0;
end

latestVersion = sscanf(latestVersion, '%d.%d.%d')';
if length(latestVersion) < 3
   latestVersion(3) = 0;
end

answer = (sign(version - latestVersion) * [10; .1; .001]) < 0;

end%upgrade
