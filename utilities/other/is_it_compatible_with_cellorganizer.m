function answer = is_it_compatible_with_cellorganizer()
% IS_IT_COMPATIBLE_WITH_CELLORGANIZER True iff the system is compatible with the CellOrganizer,
% false otherwise.

% Ivan E. Cao-Berg
%
% Copyright (C) 2008-2019 Murphy Lab
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

% March 11, 2014 I. Cao-Berg Updated method to use native isunix, ismac and
% ispc methods
%
% February 2, 2014 I. Cao-Berg Changed method to return true if the
% computer is a PC but issue a warning
%
% March 9, 2016 I. Cao-Berg Fixed logical bug that was preventing users
% from running Matlab 2015a.

if( nargin ~= 0 )
    error('CellOrganizer: Wrong number of input arguments');
else
    if ( ismac || isunix  || ispc )
        matlab_version = version;
        
        if isempty(strfind( matlab_version, '2018a' )) && ...
                isempty(strfind( matlab_version, '2018b' )) && ...
                isempty(strfind( matlab_version, '2019b' )) && ...
                isempty(strfind( matlab_version, '2019a' )) 
            warning( [ 'The current stable release of CellOrganizer has not been ' ...
                'tested with this version of Matlab. This version is ' ...
                matlab_version ] );
            answer = true;
            return
        elseif strfind( matlab_version, '2014a' )
            toolboxes = { 'Image Processing Toolbox', ...
                'Statistics Toolbox', ...
                'Curve Fitting Toolbox' };
            strfind( matlab_version, '2015a' )
            toolboxes = { 'Image Processing Toolbox', ...
                'Statistics and Machine Learning Toolbox', ...
                'Curve Fitting Toolbox' };
        else
            toolboxes = { 'Image Processing Toolbox', ...
                'Statistics and Machine Learning Toolbox', ...
                'Curve Fitting Toolbox' };
        end
        
        for index=1:1:length( toolboxes )
            toolbox = toolboxes{index};
            if ~has_toolbox( toolbox )
                warning( [ 'Unable to find toolbox: ' ...
                    toolbox '.'] );
                answer = false;
                return
            end
            answer = true;
        end
    end
    
    if ispc()
        warning( ['CellOrganizer has not been tested in any ' ...
            'Windows flavor. Some functionality might not be ' ...
            ' working.'] );
        warning('Unsupported operating system.');
        answer = false;
    end
end
end%is_it_compatible_with_cellorganizer

function answer = has_toolbox( name )
%HASTOOLBOX True iff the current system has the Matlab toolbox available,

if( nargin ~= 1 )
    error( 'CellOrganizer: Wrong number of input arguments' );
else
    answer = false;
    information = ver;
    for( i = 1:1:length( information ) )
        if( strcmpi( information(i).Name, name ) )
            answer = true;
            return;
        end
    end
end
end%has_toolbox
