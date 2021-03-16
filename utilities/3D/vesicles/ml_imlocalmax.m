function pts = ml_imlocalmax(img)
%ML_IMLOCALMAX Find local maxima in an image.
%   PTS = ML_IMLOCALMAX(IMG) returns a [point array] containing the
%   coodinates of the local maxima in the [image] IMG.
%   
%   See also

%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('Exactly 1 argument is required');
end

% Another algorithm
% statusMatrix = ones(size(img));
% statusMatrix(1,:) = 0;
% statusMatrix(:,1) = 0;
% statusMatrix(end,:) = 0;
% statusMatrix(:,end) = 0;
% 
% offset = [0 1;1 0;1 1;1 -1];
% 
% img = [img;zeros(1,size(img,2))];
% img = [img zeros(size(img,1),1)];
% 
% for i=2:size(img,1)-1
%     for j=2:size(img,2)-1
%         for k=1:size(offset,1)
%             i2 = i+offset(k,1);
%             j2 = j+offset(k,2);
%             if(img(i,j)<img(i2,j2))
%                 statusMatrix(i,j) = 0;
%             else
%                 statusMatrix(i2,j2) = 0;
%             end
%         end
%     end
% end


if size(img,3)<=1
    
    sub = img;
    subc = sub(2:end-1,2:end-1);
    subl = sub(2:end-1,1:end-2);
    subr = sub(2:end-1,3:end);
    subu = sub(1:end-2,2:end-1);
    subd = sub(3:end,2:end-1);
    subul = sub(3:end,1:end-2);
    subur = sub(1:end-2,3:end);
    subll = sub(1:end-2,1:end-2);
    sublr = sub(3:end,3:end);

    statusMatrix = (subc > subl) & (subc > subr) & ...
        (subc > subu) & (subc > subd) & ...
        (subc > subul) & (subc > subur) & (subc > subll) & (subc > sublr);
    
    if all(statusMatrix==0)
        pts = [];
    else
        [pts(:,1),pts(:,2)] = find(statusMatrix==1);
        pts = pts+1;
    end

elseif size(img,3)>1

%----------------------------------
%3D version
%   18-Jan-2007 Written by T. Peng

    [L,H,D] = size(img);
    subc = img(2:L-1,2:H-1,2:D-1);
    k = 1;  statusMatrix = true((L-2)*(H-2)*(D-2),1);
    for l=-1:1
        for h=-1:1
            for d=-1:1
                if l||h||d
                    sub = img(2+l:L-1+l,2+h:H-1+h,2+d:D-1+d);
                    compare=subc>sub;
                    statusMatrix=statusMatrix & compare(:);
                end
            end
        end
    end

    if all(statusMatrix==0)
        pts = [];
    else
        [pts(:,1),pts(:,2),pts(:,3)] = ...
            ind2sub([L-2,H-2,D-2],find(statusMatrix==1));
        pts = pts + 1;
    end
%---------------------------------

end