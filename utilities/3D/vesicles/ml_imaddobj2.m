function [img2,obj_keep] = ml_imaddobj2(img,objs,param, options)
%ML_IMADDOBJ2 Add an object into an image.
%   IMG2 = ML_IMADDOBJ2(IMG,OBJ) adds an [object] into the [image] IMG and
%   returns the new image. The position of OBJ is determined by its 
%   coordinates.
%   
%   IMG2 = ML_IMADDOBJ2(IMG,OBJ,PARAM) customizes different ways of putting
%   OBJ. PARAM can have the following fields:
%   ------------------------------------------------------------------
%        name      | data type |  description           | Default value
%   ------------------------------------------------------------------
%      pos         | [point]   | position of OBJ        |      []
%     method       | string    | how to add OBJ         |   'replace' 
%      posref      | string    | how to specify object  |   'cof'
%                  |           | position               |
%   ------------------------------------------------------------------
%
%   More details:
%       'pos' - if it is empty, the position of the oject will be
%           determined by its coordinates.
%       'method' - it has two values. 'replace' means that the the values
%           in IMG will be replaced by OBJ intensities. 'add' means that
%           the OBJ intensities will be add into IMG.
%       'posref' - it is only applicable when 'pos' is not empty. If it is
%           'cof', then 'pos' will be where the COF of OBJ is located in
%           the image IMG. If it is 'corner', then 'pos' will be where the
%           left top coner of the bounding box of OBJ is located.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU
%   25-April 2012 Added objectmethod detection to speed up calculations:DPS 
% August 7, 2012 D. Sullivan Reverted back to adding objects to the images
%                            after they were all sampled. 
% August 7, 2012 D. Sullivan Created check for whether the obj list is a
%                            single object of a cell of objects
% August 22, 2013 G. Johnson 'Objectmethod' no longer a manditory parameter
%                            field
%
% Copyright (C) 2007  Murphy Lab
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

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('method','replace','pos',[],'posref','cof', 'objectmethod', []));

%adding slice above and below to check for oob vesicles
cell2 = false(size(options.cell) + [0,0,2]);
cell2(:,:,2:end-1) = options.cell;

img2 = zeros(size(img) + [0,0,2], class(img));
img2(:,:,2:end-1) = img;
imageSize = size(img2);
%%
%added by DPS 3/18/12
if iscell(objs)
  objcell = objs;
end
obj_keep = false(size(param.pos, 1));
for i = 1:size(param.pos,1)
    if exist('objcell','var')
      obj = objcell{i};
    end
    
    %account for extra slice added to image size for oob check
    obj(:,3) = obj(:,3) + 1;

    if ~isempty(param.pos)
       switch param.posref
          case 'cof'
              offset = ml_calcobjcof(obj);
          case 'corner'
              box = ml_boundbox(obj(:,1:3));
              offset = box(1,:);
       end
      obj(:,1:3) = round(ml_addrow(obj(:,1:3),-offset+param.pos(i,:)));
    end



    idx = ml_objinimg(obj,imageSize);
    if ~isempty(idx)
        %warning('The object is out of the image range!')
        obj(idx,:) = [];
    end
    objidx = sub2ind(imageSize,obj(:,1),obj(:,2),obj(:,3));

    if ~isempty(obj)
      switch param.method
          case 'replace'
              %added disc object method if statement
              %DPS 4/25/12
              if (strcmpi(param.objectmethod,'disc'))
                  obj(:,4) = 255;
              end
                  
            if mean(cell2(objidx)) >= (1-options.oobthresh)  
                img2(objidx) = obj(:,4);
                obj_keep(i) = true;
            end
            img2 = img2 .* cell2;
          case 'add'
%             for k=1:size(img2,3)
%                 img3 = img2(:,:,k);
%                 img3(objidx) = img3(objidx)+obj(:,3);
%                 img2(:,:,k) = img3;
%             end
            %added sampling method within this DPS 4/26/12
            if isfield(param, 'objectmethod') & strcmpi(param.objectmethod,'sampled')
              
              %zero out object intensities in preperation for sampling
              obj(:,4) = 0;
  
              %get the number of samples 
              try
                  if iscell(param.numsamplescell)
                      numsamples = param.numsamplescell{i};
                  else
                      numsamples = param.numsamplescell;
                  end
              catch
                error('Number of samples not specified');
              end
              
              %generate positions of molecules
              %Note: multiple molecules can be assigned to the
              %"same" position due to resolution
              Rnumrnd = randi(length(obj),1,numsamples);
   
              %figure out counts of molecules in each position
              [counts,positions] = hist(Rnumrnd,length(obj));
              
              obj(:,4) = counts';
            end
            
            %add the object to the image
            if sum(imag(objidx))~=0
                keyboard
            end
            if sum(imag(obj(:,4)))~=0
                keyboard
            end
            img2(objidx) = img2(objidx)+obj(:,4);
            obj_keep(i) = true;
          otherwise
              error(['Unrecognized object adding method: ' param.method]);
      end
    end
end
%remove extra slices added for oob check
img2 = img2(:, :, 2:end-1);
%end 3/18/12 change
%%
