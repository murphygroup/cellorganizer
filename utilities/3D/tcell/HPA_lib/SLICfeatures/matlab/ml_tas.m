function [names, values, slfnames] = ml_tas(img,original_version)
% ML_TAS --- Calculate Threshold Adjacency Statistics
% [NAMES, VALUES, SLFNAMES] = ML_TAS(IMAGE, ORIGINALVERSION);
%
%  This algorithm was presented by Hamilton et al.
% in "Fast automated cell phenotype image classification"
% (http://www.biomedcentral.com/1471-2105/8/110)
%
%  The current implementationis an adapted version which
% is free of parameters. The thresholding is done beforehand
% (automatically), the margin around the mean of pixels to be
% included is the standard deviation of the pixel values
% and not fixed to 30, as before.
%
%  Set original_version to true to get the original version of the features.
% Also do NOT run the preprocessing code on images on which the original
% are to be calculated on.
% 
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

if nargin < 2,
    original_version = 1;
end

pixels = double(img(:));
if original_version,
    pixels(pixels <= 30) = [];
else,
    pixels(pixels == 0)=[];
end

mu = mean(pixels);

if original_version,
    Margin=30;
else,
    Margin=std(pixels);
end

%disp(['mu: ' num2str(mu)])
%disp(['margin: ' num2str(Margin_Top)])

bimg=(img > mu - Margin) .* (img < mu + Margin);

values=calculatetas(1.-bimg);
values=[values calculatetas(bimg)];

names = {'tas_0', 'tas_1', 'tas_2', 'tas_3', 'tas_4', 'tas_5', 'tas_6', 'tas_7', 'tas_8', ...
    'ntas_0', 'tas_1', 'tas_2', 'tas_3', 'tas_4', 'tas_5', 'tas_6', 'tas_7', 'tas_8' };
slfnames = {'tas_0', 'tas_1', 'tas_2', 'tas_3', 'tas_4', 'tas_5', 'tas_6', 'tas_7', 'tas_8', ...
    'ntas_0', 'tas_1', 'tas_2', 'tas_3', 'tas_4', 'tas_5', 'tas_6', 'tas_7', 'tas_8' };
end

function [values] = calculatetas(bimg),
M=[1 1 1; 1 10 1; 1 1 1];
V=conv2(double(bimg),double(M),'valid');
values=imhist(uint8(V));
values=values(1:9);
T=sum(values);
if T > 0,
    values=values'/T;
else,
    values=ones(size(values))'/9; % Return a fixed value instead of crashing...
end

end

% vim: set ts=4 sts=4 expandtab smartindent:
