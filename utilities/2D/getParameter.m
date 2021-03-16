function parameter = getParameter( snippet )
%GETPARAMETER Returns the information of the parameter from the XML
%snippet.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 29, 2007
% Last Update: March 8, 2008
%
% Copyright (C) 2008 Center for Bioimage Informatics/Murphy Lab
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

if( isa( snippet, 'char' ) || isa( snippet, 'cell' ) )
    parameter.id = getId( snippet );
    parameter.name = getName( snippet );
    parameter.constant = getConstant( snippet );
    parameter.complex = getComplex( snippet );
    
    if( strcmpi( parameter.complex, 'true' ) )
        parameter.type = [];
    else
        parameter.type = getType( snippet );
    end
    
    parameter.notes = getNotes( snippet );
    parameter.value = getValue( snippet, parameter.type );
else
    error( 'CellOrganizer: Input argument must be a string or cell array' );
end
end%getParameter

function id = getId( snippet )
if isa( snippet, 'char' )
    delimiters = findstr( snippet, 'id' );
else
    delimiters = findstr( snippet{1}, 'id' );
end

if( isempty( delimiters ) )
    id = [];
else
    delimiters = delimiters+4;
    for i=delimiters:1:length(snippet)
        if( findstr( snippet, '"' ) )
            delimiters = [delimiters i-1];
            break
        end
    end

    id = snippet(delimiters(1):delimiters(2));
end
end%getId

function name = getName( snippet )
if isa( snippet, 'char' )
    delimiters = findstr( snippet, 'name' );
else
    delimiters = findstr( snippet{1}, 'name' );
end

if( isempty( delimiters ) )
    error('CellOrganizer: In SLML Level 1 the attribute name is required for every parameter' );
else
    if( length(delimiters) > 1 )
        delimiters(1) = delimiters(1)+6;

        for i=delimiters(1):1:length(snippet)
            if( strcmpi( snippet(i), '"' ) )
                delimiters(2) = i-1;
                break
            end
        end
    else
        delimiters = delimiters+6;
        for i=delimiters:1:length(snippet)
            if( strcmpi( snippet(i), '"' ) )
                delimiters = [delimiters i-1];
                break
            end
        end
        name = snippet(delimiters(1):delimiters(2));
    end
end
end%getName

function constant = getConstant( snippet )
if isa( snippet, 'char' )
    delimiter = findstr( snippet, 'constant' );
else
    delimiters = findstr( snippet, 'constant' );
end

if( isempty( delimiter ) )
    constant = 'true';
else
    delimiter = delimiter+10;
    if( strcmpi( delimiter, 't' ) )
        constant = 'true';
    else
        constant = 'false';
    end
end
end%getConstant

function complex = getComplex( snippet )
if isa( snippet, 'char' )
    delimiter = findstr( snippet, 'complex' );
else
    delimiter = findstr( snippet{1}, 'complex' );
end

if( isempty( delimiter ) )
    complex = 'false';
else
    delimiter = delimiter+9;
    if( strcmpi( snippet(delimiter), 't' ) )
        complex = 'true';
    else
        complex = 'false';
    end
end
end%getComplex

function type = getType( snippet )
if isa( snippet, 'char' )
    delimiter = findstr( snippet, 'type' );
else
    delimiter = findstr( snippet{1}, 'type' );
end

if( isempty( delimiter ) )
    type = 'double';
else
    delimiter = delimiter+6;
    if( strcmpi( snippet(delimiter), 'd' ) )
        type = 'double';
    elseif( strcmpi( snippet(delimiter), 'i' ) )
        type = 'integer';
    elseif( strcmpi( snippet(delimiter), 's' ) )
        type = 'string';
    else
        type = 'boolean';
    end
end
end%getType

function notes = getNotes( snippet )
if isa( snippet, 'char' )
    delimiters = findstr( snippet, 'notes' );
else
    delimiters = findstr( snippet{1}, 'notes' );
end

if( isempty( delimiters ) )
    notes = [];
else
    delimiters = delimiters+7;
    for i=delimiters:1:length(snippet)
        if( findstr( snippet, '"' ) )
            delimiters = [delimiters i-1];
            break
        end
    end
    notes = snippet(delimiters(1):delimiters(2));
end
end%getsNotes

function value = getValue( snippet, type )
if isa( snippet, 'char' )
    delimiters = findstr( snippet, 'value' );
    if( isempty( delimiters ) )
        value = [];
    else
        delimiters = delimiters+7;
        for i=delimiters:1:length(snippet)
            if( findstr( snippet, '"' ) )
                delimiters = [delimiters i-1];
                break
            end
        end

        value = snippet(delimiters(1):delimiters(2));
        
        switch type
            case{'double', 'integer'}
                value = str2double(value);
            case 'boolean'
                if( strcmpi( value, 'true' ) )
                    value = true;
                else
                    value = false;
                end
            otherwise
                value = char(value);
        end
    end
else
    value = mathml2array( snippet(2:end-1) );
end
end%getValue



