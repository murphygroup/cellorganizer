function setup( force )
%SETUP Helper method that loads CellOrganizer into workspace if toolbox is
%compatible. If force flag is set to true, then it will load CellOrganizer
%even if system is not compatible.

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2008-2020 Murphy Lab
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

%icaoberg 2/3/2015
if ~exist( 'force', 'var' )
    force = false;
end

addpath( pwd );
global CELLORGANIZER_ROOT_PATH
CELLORGANIZER_ROOT_PATH = pwd;

%icaoberg 10/2/2015
disp('Adding appropiate folders to path.');
addpath(genpath( [pwd filesep 'utilities']));
addpath(genpath( [pwd filesep 'models']));
addpath(genpath( [pwd filesep 'demos']));
addpath(genpath( [pwd filesep 'applications']));
addpath([pwd filesep 'images']);

javaaddpath( [ pwd filesep ...
    'utilities' filesep 'bfmatlab' filesep 'bioformats_package.jar'] );

%enable logging
loci.common.DebugTools.enableLogging('INFO');

%icaoberg 2/3/2015
disp('Checking if your system and Matlab version is compatible with CellOrganizer.');
if ~is_it_compatible_with_cellorganizer()
    if force
        warning( ['The current system is not ' ...
            'compatible with CellOrganizer. However, CellOrganizer was loaded because force flag was set to true.'] );
    else
        rmpath(genpath( [pwd filesep 'utilities']));
        rmpath(genpath( [pwd filesep 'models']));
        rmpath(genpath( [pwd filesep 'demos']));
        error( ['The current system is not ' ...
            'compatible with CellOrganizer. Please refer to the official documentation for more information.'] );
    end
end

version = '2.9.3';
versionURL = 'http://murphylab.web.cmu.edu/software/CellOrganizer/version';

try
    fprintf( 1, '%s', 'Checking for updates. ' );
    latestVersion = urlread( versionURL );
    if latestVersion(end) == 10
        latestVersion = latestVersion(1:end-1);
    end
    
    if upgrade( version, latestVersion )
        disp(['A newer version (CellOrganizer v' latestVersion ...
            ') is available online. Please visit www.cellorganizer.org to download the latest stable release.']);
    else
        disp(['CellOrganizer version ' num2str(version) ' is the latest stable release.']);
    end
catch
    disp('Unable to connect to server. Please try again later.');
end
end%setup

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
