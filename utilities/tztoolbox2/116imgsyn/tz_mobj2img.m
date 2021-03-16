function img=tz_mobj2img(combobjects,param)
%TZ_MOBJ2IMG Align multiple objects in an image without overlaping.
%   IMG = TZ_MOBJ2IMG(COMBOBJECTS) returns an image with multiple
%   objects in it. The objects in the cell array of [object]s COMBOBJECTS
%   are aligned in the image with their originial positions.
%
%   IMG = TZ_MOBJ2IMG(COMBOBJECTS,PARAM) lets users specify how to align
%   objects. PARAM could have the following fields:
%       'method' - how to put the objects into the image. It has the
%           following values:
%           'org' - put objects into their original positions (default).
%           'mat' - align the objects like a matrix.
%       'imgsize' - the size of the returned image
%       'isscaling' - scaling or not. only for 'method' is 'mat'. This is
%           to be compatible with old files.
%           
%       
%   UNDER CONSTRUCTION       

%   ??-???-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   15-May-2006 Modified T. Zhao
%       - add more options of align objects
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('method','org','imgsize',[512 512]));

switch(param.method)
    case 'org'
        img = tz_objs2img(combobjects,param.imgsize);
    case 'mat'
        if isscaling==1
            orgsize=[512 512];
            img=zeros(orgsize);

            curplace=[1,1];
            pro=curplace(2);

            for i=1:length(combobjects)
                subimg=tz_obj2img(combobjects{i},[]);
                subimg=subimg*255/max(subimg(:));
                subsize=size(subimg);
                if subsize(1)+curplace(1)-1<=orgsize(1)
                    if subsize(2)+curplace(2)-1<=orgsize(2)
                        img(curplace(1):(curplace(1)+subsize(1)-1), ...
                            curplace(2):(curplace(2)+subsize(2)-1))=subimg;
                        curplace(1)=curplace(1)+subsize(1);
                        pro=max(pro,curplace(2)+subsize(2));
                    end
                else
                    if subsize(2)+pro-1<=orgsize(2)
                        curplace(2)=pro;
                        curplace(1)=1;
                        img(curplace(1):(curplace(1)+subsize(1)-1), ...
                            curplace(2):(curplace(2)+subsize(2)-1))=subimg;
                        curplace(1)=curplace(1)+subsize(1);
                        pro=max(pro,curplace(2)+subsize(2));
                    end
                end
            end

            img=imresize(img,imgsize,'bilinear');
        else
            for i=1:length(combobjects)
                [box,s(i,:)] = tz_boundbox(combobjects{i});               
            end
            maxs = max(s,[],1)+param.space;
            nrow = floor(param.imgsize(1)/maxs(1));
            ncol = floor(param.imgsize(1)/maxs(2));
            k = 1;
            for i=1:nrow
                for j=1:ncol
                    
                end
            end
        end
    case 'syn'
        img=zeros(param.imgsize);
        objsizes = tz_objsizes(combobjects);
        [objsizes,idx] = sort(objsizes,'descend');
        combobjects = combobjects(idx);
        dists = param.dists(idx);
        for i=1:length(combobjects)
            overlap = 1;
            maxtrial = 20;
            trial = 0;
            while overlap>0
                theta = rand(1)*360;
                newobj = tz_rotateobj(combobjects{i},theta);
                [img,overlap] = tz_imaddobj(newobj,img,dists(i),[], ...
                    {2,'random','original'});
                trial = trial+1;
                if trial>maxtrial
                    [img,overlap] = tz_imaddobj(newobj,img, ...
                        dists(i),[],{2,'optimal','original'});
                    if overlap>0
                        warning(['An object with size '  ...
                            num2str(size(newobj,1)) ...
                            ' was failed to put in the image']);
                        break;
                    end
                end
            end
        end
end

