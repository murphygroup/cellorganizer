function threshold = ml_rcthreshold(img)  
% ML_RCTHRESHOLD - Find a threshold for IMAGE 
% ML_RCTHRESHOLD(IMAGE) - input is output from nihscale(image)
% output is threshold value using same algorithm as nih image
% Ridler and Calvard (IEEE tans. on Systems, man, cybernetics
% SMC-8 no 8 aug 1978)
% Trussel (same journal) SMC-9 no 5 may 1979
%
% Basically, the algorithm moves the proposed threshold from the 
% minimum non-zero bin of the histogram to the maximum non-zero bin.
% While doing so, it calculates the 'center of mass' separately
% for that part of the histogram below the proposed threshold and
% that part above.  While the proposed threshold is less than 
% the mean of these two centers of mass, the algorithm continues.
% When this condition no longer holds, the proposed threshold is
% rounded to give a final result.
%
% If <= one intensity value is nonzero output threshold will be 1/2
% and error message will be printed to screen

% Copyright (C) 2006-2012 Murphy Lab
% Carnegie Mellon University
%
% 20 Aug 98 - M.V. Boland.  Modified from version written by W. Dirks.
%
% Sometime in Year 2000, MV and RFM noticed that the algorithm
% had not been implemented properly. Corrected by MV and RFM.
%
% Modified by T. Zhao to fix bugs
% July 2, 2012 Gregory Johnson modified to accept images of arbitrary bit depth by converting them do 16-bit
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

converted = false;
if strcmp(class(img),'uint8')
    %disp ('Image is uint8, using 256 bin histogram');
    bins = 256;
    
elseif strcmp(class(img),'uint16')
    %disp ('Image is uint16, using 65536 bin histogram');
    bins = 65536;
else
    %disp ('Can''t determine bit depth. Converting to 16-bit')
    bins = 65536 ;
    
    img = double(img);
    
    
    
    imgMin = min(img(:));
    imgMax = max(img(:));
    
    imgDiff = imgMax - imgMin;
    
    
    
    img = uint16((img- imgMin)/imgDiff * bins);

    converted = true;
end

imgsize=size(img);

if length(imgsize)>2
    img=reshape(img,imgsize(1),prod(imgsize(2:end)));
end

histo = imhist(img(:),bins);

hsum = sum(histo);


histo(1) = [];
histo(bins-1) = 0;


whitepos=find(histo>0);

MinIndex=1;
MaxIndex=1;
if ~isempty(whitepos) 
    MinIndex=whitepos(1);
    MaxIndex=whitepos(end);
end

if (MinIndex >= MaxIndex)
	threshold = MinIndex;
	warning('There is <= 1 nonzero pixel intensity');
	return;
end;

MovingIndex = MinIndex+1;

Prev = MinIndex;
while (Prev ~= MovingIndex)
    if (MovingIndex >= MaxIndex)
        break;
    end;
    
    sum1 = sum([MinIndex:MovingIndex-1]' .* ...
        histo(MinIndex:MovingIndex-1)) ;
    sum2 = sum(histo(MinIndex:MovingIndex-1)) ;
    sum3 = sum([MovingIndex:MaxIndex]' .* ...
        histo(MovingIndex:MaxIndex)) ;
    sum4 = sum(histo(MovingIndex:MaxIndex)) ;
    
    Prev = MovingIndex;
    loweravg=sum1/sum2;
    higheravg=sum3/sum4;
    
    MovingIndex = ceil((loweravg + higheravg)/2) ;   
end

threshold = MovingIndex;

if converted
    threshold = (threshold/bins * imgDiff) + imgMin;
end
