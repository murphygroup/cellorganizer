function threshold = ml_threshold(image)  
% ML_THRESHOLD - Find a threshold for IMAGE 
% ML_THRESHOLD(IMAGE) - input is output from nihscale(image)
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
%

% Copyright (C) 2006  Murphy Lab
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

% 20 Aug 98 - M.V. Boland.  Modified from version written by W. Dirks.
%

% $Id: ml_threshold.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

%
% Generate a 256 bin histogram and 'flip' it.
%

if length(size(image))>2
    image=reshape(image,size(image,1),prod(size(image))/size(image,1));
end

histo = flipdim(imhist(image,256),1) ;

histo(1) = 0;
histo(256) = 0;
error = 0; 

MinIndex = 1;
MaxIndex = 256;

while ( (histo(MinIndex) == 0) & (MinIndex < 256) )
	MinIndex = MinIndex + 1;
end;

while ((histo(MaxIndex) == 0) & (MaxIndex > 1))
	MaxIndex = MaxIndex - 1;
end;

if (MinIndex >= MaxIndex)
	Result = 128;
	disp('There is <= 1 nonzero pixel intensity');
	error = 1;
end;

MovingIndex = MinIndex;

if (error == 0)
	Result = MaxIndex;
	while ((MovingIndex+1)<=Result)
		if (MovingIndex > (MaxIndex-1))
			break;
		end;

		sum1 = sum([MinIndex:MovingIndex]' .* ...
                            histo(MinIndex:MovingIndex)) ;
		sum2 = sum(histo(MinIndex:MovingIndex)) ;
		sum3 = sum([MovingIndex+1:MaxIndex]' .* ...
                            histo(MovingIndex+1:MaxIndex)) ;
		sum4 = sum(histo(MovingIndex+1:MaxIndex)) ;

		Result = (sum1/sum2 + sum3/sum4)/2 ;
                         
		MovingIndex = MovingIndex+1;
	end;

end;

level = round(Result);
% value = (256 - level);
level2 = 1 - level/255;

threshold = level2;
