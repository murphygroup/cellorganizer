function [allfiles, patternnum] = ml_ls( patterns )
% FILELISTING = ML_LS( PATTERN ) Returns file listing of wildcard PATTERN (e.g. '*.jpg').
% FILELISTING is a cell array of strings of filenames without path.
% the listing is sorted by name alphabetically
%

% 8/11/13 G. Johnson - modified such that pattern can be a string, or cell
%                      array of strings, in case of the need for multiple
%                      source directories
%
% Copyright (C) 2006-2013 Murphy Lab
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

if ~iscell(patterns)
    patterns = {patterns};
end

allfiles = cell(length(patterns),1);

%hidden option to the user
option = 1;

if option == 1
    for c = 1:length(patterns)
        try
            files = {};
            [~, list] = system(['ls -1 ' patterns{c}]);
            rows = strread( list, '%s', 'delimiter', '\n');
            for i=1:1:length(rows)
                temp = {};
                %xruan 06/04/2015
                %add the judge of empty directory with system command
                if  isempty(strfind(list, 'No such file or directory'))
                    temp = strread(rows{i},'%s','delimiter','\t');
                end
                for j=1:1:length(temp)
                    if ~isempty(temp{j})
                        files{length(files)+1} = temp{j};
                    end
                end
            end
        catch
            files = {};
        end
        
        %D. Sullivan 6/20/13 - added natural sorting
        files = sort_nat( files );
        allfiles{c} = files;
        patternnum{c} = ones(1,length(files))*c;
    end
    
    allfiles = [allfiles{:}];
    patternnum = [patternnum{:}];
else
    for c = 1:length(patterns)
        try
            files = {};
            [~, list] = system(['ls -1 ' patterns{c}]);
            rows = strread( list, '%s', 'delimiter', '\n');
            for i=1:1:length(rows)
                temp = {};
                %xruan 06/04/2015
                %add the judge of empty directory with system command
                if  isempty(strfind(list, 'No such file or directory'))
                    temp = strread( rows{i},'%s','delimiter','\t');
                end
                for j=1:1:length(temp)
                    if ~isempty(temp{j})
                        files{length(files)+1} = temp{j};
                    end
                end
            end
        catch
            files = {};
        end
        
        %D. Sullivan 6/20/13 - added natural sorting
        files = sort_nat( files );
        allfiles{c} = files;
        patternnum{c} = ones(1,length(files))*c;
    end
    
    allfiles = [allfiles{:}];
    patternnum = [patternnum{:}];
end
end
