function img2 = tz_imtranslate(img,offset,option)
%TZ_IMTRANSLATE Translate an image.
%   IMG2 = TZ_IMTRANSLATE(IMG,OFFSET) returns an image that is the
%   translation of IMG with offset OFFSET, which is a 1x2 vector.
%   
%   IMG2 = TZ_IMTRANSLATE(IMG,OFFSET,OPTION) also let the user specifies
%   the post processing of translation:
%       'blank' - leave the margin blank
%       'circular' - circular shift
%
%   See also ML_TRANSLATE

%   19-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('2 or 3 arguments are required')
end

if nargin < 3
    option = 'blank';
end

switch option
case 'blank'
    img2 = shift(img,offset);
    if offset(1)~=0
        if offset(1)>0
            blankRows = 1:offset(1);
        else
            blankRows = (size(img,1)+offset(1)+1):size(img,1);
        end
        blankRows(blankRows<0) = [];
        if ~isempty(blankRows)
            img2(blankRows,:) = 0;
        end
    end
    
    if offset(2)~=0
        if offset(2)>0
            blankColumns = 1:offset(2);
        else
            blankColumns = (size(img,2)+offset(2)+1):size(img,2);
        end
        blankColumns(blankColumns<0) = [];
        if ~isempty(blankColumns)
            img2(:,blankColumns) = 0;
        end
    end  
case 'circular'
    img2 = shift(img,offset);
case 'bilinear'
    [xgrid,ygrid] = meshgrid(1:size(img,1),1:size(img,2));
    xi = xgrid-offset(1);
    yi = ygrid-offset(2);
    
    img2 = ...
        interp2(xgrid,ygrid,double(img'),xi,yi,'bilinear');
    img2 = img2';
    img2(isnan(img2)) = 0;
case 'otherwise'
    error('Invalid option');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x2 = shift(x, offset)

dims = size(x);

offset = mod(-offset,dims);

x2 = [ x(offset(1)+1:dims(1), offset(2)+1:dims(2)),  ...
    x(offset(1)+1:dims(1), 1:offset(2)); ...
    x(1:offset(1), offset(2)+1:dims(2)), ...
    x(1:offset(1), 1:offset(2)) ];
        