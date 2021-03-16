function stack = ml_tiffread( tifffile, channels, xrange, yrange, zrange, tmpdir)

% IMAGE = ML_TIFFREAD( TIFFFILE, CHANNELS, XRANGE, YRANGE, ZRANGE)
%
% Reads a 3D image from a TIFF file.
% TIFFFILE is the handle of the tiff file (returned by ml_tiffopen).
% CHANNELS is a vector of numbers (e.g. [2 3 1]) indicating which
% channels to load in which order.
% XRANGE, YRANGE and ZRANGE specify the region to be loaded.
% e.g. im = readtiff(tiff_f,[1 3],[600:900],[800:950],[5:24]);
% An empty matrix [] specifies full range (omission of arguments
% has the same effect).

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

%%% PROCESS THE ARGUMENTS
if( nargin > 6) error('ML_TIFFREAD takes at most 6 arguments'); end
if( nargin < 6) tmpdir = tempdir; end
if( nargin < 1) error('ML_TIFFREAD takes at least one argument - the TIFF file handle'); end
if( nargin < 5) zrange = []; end
if( nargin < 4) yrange = []; end
if( nargin < 3) xrange = []; end
if( nargin < 2) channels = []; end

% If TIFFFILE is a string, then use mv_tiffopen first
if( ischar(tifffile))
    switch tifffile(end-2:end)
     case '.gz'
         is_compressed = 1;
         uncompress_cmd = 'gunzip -c';
     case 'bz2'
         is_compressed = 1;
         uncompress_cmd = 'bunzip2 -c';
     otherwise
         is_compressed = 0;
    end
    if( is_compressed
        [f,tempfilename] = ml_fopentemp(tmpdir);
        tempfullpath = [tmpdir tempfilename];
        fclose(f);
        unix([uncompress_cmd ' ' tifffile ' > ' tempfullpath]);
        tifffile = tempfullpath;
    end
    tifffile = ml_tiffopen( tifffile);
end

% Deal with unspecified arguments
if( isempty(zrange))
    no_z_range = 1;
else
    no_z_range = 0;
end;
if( isempty(xrange) & isempty(yrange))
    no_xy_range = 1;
else
    if( isempty(xrange) | isempty(yrange))
        error(['XRANGE and YRANGE can only be specified together. '...
	       'You cannot specify one and then not specify the' ...
	       ' other']);
    else
        no_xy_range = 0;
    end
end;
n_chans = length(tifffile.ChannelInfo);
if( isempty(channels))
    channels = [1 : n_chans];
else
    if( is_out_of_bounds( channels, [1 n_chans]))
      error('Invalid Channel No.-s specified in CHANNELS');
    end
end

% Open the file
f = fopen(tifffile.Filename,'rb',tifffile.Endianity);
% Read the images from each channel
NoOfChannels = length(channels);
%if( tifffile.ChannelInfo{channels(1)}.PlanarConfig == 1) % interleaved
%  NoOfReadsPer
%else
%end
for i = 1 : NoOfChannels
  c = channels(i);
  BytesPerSample = ceil(tifffile.ChannelInfo(c).BitsPerSample/8);
  if( length( BytesPerSample) > 1)
    if( ~(BytesPerSample(1) == BytesPerSample(2:end)))
      error('BitsPerSample has different values for different channels');
    end
    BytesPerSample = BytesPerSample(1);
  end
  switch BytesPerSample
   case 1
      Voxel = uint8(0);
      VoxelType = 'uint8';
   case 2
      Voxel = uint16(0);
      VoxelType = 'uint16';
   case {3,4}
      warning('Pixels per Sample is greater than 16. Not sure if this works');
      Voxel = uint32(0);
      VoxelType = 'uint32';
   otherwise
      error('We do not handle BitsPerPixel values higher than 32.');
  end

  fprintf( 1, 'Loading channel %i: ', c);
  Height = tifffile.ChannelInfo(c).Height;
  Width = tifffile.ChannelInfo(c).Width;
  Depth = tifffile.ChannelInfo(c).Depth;
  % Subregion selection
  if( no_z_range);
      zrange = [1 : Depth]; 
      FinalDepth = Depth;
  else
      if( is_out_of_bounds( zrange, [1 Depth]))
	error('ZRANGE is out of bounds');
      end
      FinalDepth = length(zrange);
  end;
  if( no_xy_range)
      FinalWidth = Width;
      FinalHeight = Height;
  else
      if( is_out_of_bounds( xrange, [1 Width]) | ...
	  is_out_of_bounds( yrange, [1 Height]))
	error('XRANGE or YRANGE is out of bounds');
      end
      FinalWidth = length(xrange);
      FinalHeight = length(yrange);
  end;
  % Pre-create matrix to receive loaded image
  stack{i} = repmat(Voxel,[FinalHeight FinalWidth FinalDepth]);
  % do the loading
  for z = 1 : length(zrange);
    d = zrange(z);
    imageplane = ReadPixelData( f, ...
				  tifffile.ChannelInfo(c).StripOffsets(d,:), ...
				  tifffile.ChannelInfo(c).StripByteCounts(d,:), ...
				  Width, ...
				  Height, ...
				  BytesPerSample, ...
				  tifffile.ChannelInfo(c).PlanarConfig, ...
                                  tifffile.ChannelInfo(1).ImageType, ...
				  Voxel, ...
				  VoxelType);
    if( no_xy_range)
    else %%% If user requested a particular region of the image, obtain the appropriate region
        for r = 1:length(imageplane); imageplane{r} = imageplane{r}(yrange,xrange); end
    end
    switch( tifffile.ChannelInfo(1).ImageType)
      case 'RGB', %%% If image is RGB then all three channels will have been read from the same IFD
	stack{1}(:,:,z) = imageplane(:,:,1);
	stack{2}(:,:,z) = imageplane(:,:,2);
	stack{3}(:,:,z) = imageplane(:,:,3);
      otherwise,
        stack{i}(:,:,z) = imageplane(:,:,1);
    end

    fprintf( 1, '%i ', d);
  end
  fprintf( 1, '\n');
end
fclose( f);

% If there is only a single channel then do not use cell array
if( length(stack) == 1), stack = stack{1}; end;

% Remove the temporary file if one was created
if( exist('tempfullpath','var'))
    unix(['rm -f ' tempfullpath]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = ReadPixelData( f, offsets, bytecounts, ImageWidth, ImageLength, ...
			      BytesPerSample, PlanarConfig, ImageType, Voxel, VoxelType)

 %RowsPerStrip = bytecounts(1)/(BytesPerSample*ImageWidth);
 %img = repmat(Voxel,[ImageWidth ImageLength 3]); %%% 3 is for RGB
 img = [];
 for i = 1:length(offsets)
    fseek( f, offsets(i), -1);
    SamplesPerStrip = bytecounts(i)/BytesPerSample;
    switch BytesPerSample
     case 1, Strip = uint8(fread( f, SamplesPerStrip, VoxelType));
     case 2, Strip = uint16(fread( f, SamplesPerStrip, VoxelType));
     case {3,4}, Strip = uint32(fread( f, SamplesPerStrip, VoxelType));
     otherwise, error(['Something wrong with BytesPerSample (should' ...
		       ' be 1, 2, 3 or 4)']);
    end
    img = [img; Strip];
 end
 switch ImageType
  case 'RGB',
   switch PlanarConfig
    case 1,
     img = reshape( img, [3 ImageWidth ImageLength]);
     img = shiftdim( img, 1);
    case 2,
     img = reshape( img, [ImageWidth ImageLength 3]);
    otherwise, error('PlanarConfig should be either 1 or 2');
   end
  otherwise,
   img = reshape( img, [ImageWidth ImageLength]);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = is_out_of_bounds( region, limits)
  
    mismatch = [find(region > limits(2)) find(region < limits(1))];
    if( isempty(mismatch))
      result = 0;
    else
      result = 1;
    end
