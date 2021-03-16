function answer = mat2SBMLArrayData( array, oneLine )
% MAT2SBMLARRAYDATA Simple helper function that takes an array and returns
% a string representation of spatial:arrayData that is compatible with
% libSBML

% Copyright (C) 2016-2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
% Copyright (C) 2018 Taraz Buck
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

if nargin < 2
    oneLine = true;
end
if oneLine
    array = array(:);
    array = reshape( array, [], 1 );
end
answer = mat2str( array );
if oneLine
    answer = strrep(answer, ';',' ');
else
    answer = strrep(answer, ';',sprintf('\n'));
end
answer = answer(2:end-1);
