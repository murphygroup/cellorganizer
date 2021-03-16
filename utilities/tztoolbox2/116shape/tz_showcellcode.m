function h = tz_showcellcode(cellcode,option,varargin)
%TZ_SHOWCELLCODE Show a coded cell.
%   TZ_SHOWCELLCODE(CELLCODE,OPTION) shows the cell code CELLCODE according
%   to the option OPTION:
%       'cntr' - contour
%       'hitpts' - radiate points
%       'hitxy' - relationship between nucleus shape and cell shape
%
%   TZ_SHOWCELLCODE(CELLCODE,'app',OPTION2,ORDER) shows the deformation
%   shape. OPTION2 and ORDER are the last two arguments for
%   TZ_ESTRANSFNUCELL.
%
%   TZ_SHOWCELLCODE(CELLCODE,'img',IMGSIZE) shows the cell in an image with
%   size IMGSIZE.
%
%   H = TZ_SHOWCELLCODE(...) returns the handle of the shown figure or an
%   image if OPTION is 'img'.
%
%   See also

%   06-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

    
if nargout>0 & ~strcmp(option,'img') & ~strcmp(option,'hitimg')
    h = figure;
end

switch option
    case 'cntr'
        ml_showpts_2d(cellcode.cellcontour,'ln',1);
        hold on
        ml_showpts_2d(cellcode.nuccontour,'ln',1);
        hold off
    case 'hitpts'
        ml_showpts_2d(cellcode.nucellhitpts,'ln',1);
        hold on
        ml_showpts_2d(cellcode.nuchitpts,'ln',1);
        hold off
    case 'hitxy'
        plot(cellcode.nuchitpts(:,1),cellcode.cellhitpts(:,1));
        hold on
        plot(cellcode.nuchitpts(:,2),cellcode.cellhitpts(:,2),'r-');
        hold off
    case 'app'
        A = tz_esttransfnucell(cellcode,varargin{:});
        pts = tz_deformpts_2d( ...
            ml_addrow(cellcode.nuchitpts,-cellcode.nuccenter),A);
        tz_showcellcode(cellcode,'hitpts');
        hold on
        ml_showpts_2d(ml_addrow(pts,cellcode.nuccenter),'ln',-1);
        hold off
    case 'img'
        if exist('h','var')
            close(h);
        end
        
        imgsize = varargin{1};
        
        [cellBoundbox,cellBoxsize] = ml_boundbox(cellcode.nuccontour);
        offset = round((imgsize-cellBoxsize)/2)-cellBoundbox(1,:);
        
        img = ml_obj2img( ...
            ml_addrow([cellcode.cellcontour;cellcode.nuccontour],offset), ...
            imgsize);
        se=strel('disk',4,4);
        img = imdilate(img,se);
        imshow(img,[]);
        h = img;
    case 'hitimg'
        cellcode.cellcontour = ml_showpts_2d(cellcode.nucellhitpts,'ln',1);
        cellcode.nuccontour = ml_showpts_2d(cellcode.nuchitpts,'ln',1);
        
        h = tz_showcellcode(cellcode,'img',varargin{:});
    otherwise
        error(['Unrecognized option: ' option]);
end

% if isfield(cellcode,'cellmangle')
%     title(['Major angle of the cell: ' num2str(cellcode.cellmangle)]);
% end

