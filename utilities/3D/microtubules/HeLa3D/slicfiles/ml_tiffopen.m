function info = ml_tiffopen( filename)

% TIFFFILE_HANDLE = ML_TIFFOPEN( FILENAME)
%
% If used without a return argument then displays information
% about the TIFF file, such as image size and type, number of
% channels etc.
%
% TIFFFILE_HANDLE is a structure containing detailed information
% about the file. Note that there is no need to close the file
% opened with this function.
% This info is intended to be sufficient for a program
% to read image data from the tiff file.
% * Height - 1st dimension (y)
% * Width - 2nd dimension (x)
% * Depth - 3rd dimension (z)
% * BytesPerSample
% * PlanarConfig - Interleaved or "chuncky" (1), One channel after another (2)
% * NoOfChannels
% * ImageType (e.g. 'Grayscale','RGB','Indexed Color')
% * IFDOffsets - matrix of IFD offsets. Each row is a
%                 vector of IFD offsets for one channel.
%                 The number of rows is equal to the
%                 number of channels.
% * StripOffsets - cell array with one cell per channel. Each cell
%                  contains a matrix with one row of offset values
%                  per IFD.
% * StripByteCounts - cell array with the same structure as StripOffsets
%                     (tells you the number of bytes in each strip).
  
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

% This Program was based directly on the TIFF 6.0 specifications of
% June 3, 1992.
%
% Copyright (c) 2000, 2001 Meel Velliste, Carnegie Mellon University

%%% CHECK ARGUMENTS
if( nargin < 1) error('This function needs a filename as an argument'); end;
if( ~isstr(filename)) error('This function needs a filename as an argument'); end;

%%% Get this file's full path
if( filename(1) == '/')
else
    filename = [pwd '/' filename];
end

%%% READ HEADER %%%
f = fopen( filename, 'rb');
if( f == -1) error(['Could not open file: ' filename]); end;
e = char(fread( f, 2, 'uchar'))';
switch e
 case 'II', Endianity = 'ieee-le'; %%% LITTLE ENDIAN
 case 'MM', Endianity = 'ieee-be'; %%% BIG ENDIAN
otherwise, error( ['File ' filename ' does not conform to TIFF format (error_code=1)']); 
end;
fclose( f);

% Reopen the file with the correct endianity
f = fopen( filename, 'rb', Endianity);
fseek( f, 2, -1);
% Read the magic number
MagicTIFF_ID = fread( f, 1, 'uint16');
if( MagicTIFF_ID ~= 42)
	error( ['File ' filename ' has the wrong Magic ID (should be 42)']); 
end
% Read the offset of the first IFD.
IFD_Offset = fread( f, 1, 'uint32');

% Read info from all IFD-s
IFD_count = 0;
IFDOffsets = [];
ColorMapPositions = [1];
while( IFD_Offset ~= 0)
  IFDOffsets = [IFDOffsets; IFD_Offset];
  % Read next IFD
  [IFD, ColorMap] = ReadIFD( f, IFD_Offset);
  IFD_count = IFD_count + 1;
  IFDs{IFD_count} = IFD;
  % Test if ColorMap changed from previous to this IFD
  if( IFD_count == 1) PrevColorMap = ColorMap; end
  if( ColorMap == PrevColorMap)
  else
    ColorMapPositions = [ColorMapPositions; IFD_count];
  end
  PrevColorMap = ColorMap;
  % Go on to next IFD
  IFD_Offset = IFD.NextIFDOffset;
end
% Close the file
fclose( f);

% Check whether image types of all IFDs are the same
FirstImageType = IFDs{1}.ImageType;
for i = 2 : IFD_count
  if( length(IFDs{i}.ImageType) == length(FirstImageType))
    if( IFDs{i}.ImageType == FirstImageType), break; end
  end
  warning('Not all IFDs have the same image type');
end

% Find out the number of Channels (and their locations)
switch( IFDs{1}.ImageType)
 case 'RGB',
  % If image is RGB then all three channels are in each IFD
  % rather than spread over several IFDs - therefore for file
  % reading purposes we consider it as one channel
  ChannelPositions = 1;
  NoOfChannels = 1;
 otherwise,
  % Figure out Channel numbering from info in IFDs
  % Check to see if channel numbering is included in the DocumentName
  ChannelNumbers = [];
  for i = 1 : IFD_count
    ChNo = find_channel_no( IFDs{i}.DocumentName);
    if( ChNo >= 0)
      ChannelNumbers(i) = ChNo;
    end
  end
  % If channel numbers okay in document name keep them as they are
  if( length( ChannelNumbers) == IFD_count)
    % otherwise replace by dummy
  else ChannelNumbers = zeros(1,IFD_count);
  end
  
  % Also use other attributes to find the different channels
  ChannelPositions = ColorMapPositions;
  for i = 2 : IFD_count
    if( start_next_channel( IFDs{i}, IFDs{i-1}, ...
	                    ChannelNumbers(i), ChannelNumbers(i-1)))
      if( find(i == ColorMapPositions))
      else
	ChannelPositions = [ChannelPositions; i];
      end
    end
  end
  ChannelPositions = sort(ChannelPositions);
  NoOfChannels = length( ChannelPositions);
end

% Construct info structures for each channel
ChannelPositions = [ChannelPositions; IFD_count+1];
for c = 1 : NoOfChannels
  ChanBegin = ChannelPositions(c);
  ChanEnd = ChannelPositions(c+1)-1;
  ImageDepth = ChanEnd - ChanBegin + 1;
  % Construct Strip lists for this channel
  NoOfStrips = length(IFDs{ChanBegin}.StripOffsets);
  StripOffsets = zeros(ImageDepth,NoOfStrips);
  StripByteCounts = zeros(ImageDepth,NoOfStrips);
  for i = 1 : ImageDepth
      offsets = IFDs{ChanBegin+i-1}.StripOffsets;
      bytecounts = IFDs{ChanBegin+i-1}.StripByteCounts;
      if( length(offsets) == NoOfStrips & length(bytecounts) == NoOfStrips)
	StripOffsets(i,:) = offsets;
	StripByteCounts(i,:) = bytecounts;
      else
	error(['Images of the same channel have different numbers' ...
	       ' of strips. This is bad and we do not deal with it']);
      end
  end
  chaninfo(c) = struct( 'Height', IFDs{ChanBegin}.ImageLength, ...
		    'Width', IFDs{ChanBegin}.ImageWidth, ...
		    'Depth', ImageDepth, ...
		    'ImageType', IFDs{ChanBegin}.ImageType, ...
		    'BitsPerSample', IFDs{ChanBegin}.BitsPerSample, ...
                    'PlanarConfig', IFDs{ChanBegin}.PlanarConfig, ...
		    'IFDOffsets', IFDOffsets(ChanBegin:ChanEnd), ...
		    'StripOffsets', StripOffsets, ...
		    'StripByteCounts', StripByteCounts ...
		    );
end

info = struct('Filename',filename, ...
	      'Endianity',Endianity, ...
	      'ChannelInfo',chaninfo ...
	      );

% Print out info if necessary
%fprintf( 1, '\n');
%message = sprintf( 'Total of %d slices\n', frame_no);
%fprintf( 1, message);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decision = start_next_channel( IFD, PrevIFD, ChNo, PrevChNo)
  % When certain attributes change from one IFD to the next
  % (such as CHANNEL no. in DocumentName, or ColorMap)
  % it means a new channel starts at that point

  if( IFD.ImageWidth == PrevIFD.ImageWidth) else decision = 1; return, end;
  if( IFD.ImageLength == PrevIFD.ImageLength) else decision = 1; return, end;
  if( IFD.ImageType == PrevIFD.ImageType) else decision = 1; return, end;
  if( IFD.BitsPerSample == PrevIFD.BitsPerSample) else decision = 1; return, end;
  if( IFD.SamplesPerPixel == PrevIFD.SamplesPerPixel ) else decision = 1; return, end;
  if( ChNo == PrevChNo ) else decision = 1; return, end;
  %if( IFD. == PrevIFD. ) else decision = 1; return, end;
  decision = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IFD, ColorMap] = ReadIFD( fileid, offset)

ImageWidth = [];
ImageLength = [];
BitsPerSample = [];
SamplesPerPixel = 1;
Compression = [];
PhotometricInterpretation = [];
StripOffsets = [];
StripByteCounts = [];
ColorMapOffset = 0;
ColorMap = 0;
DocumentName = [];

fseek( fileid, offset, -1);
EntryCount = fread( fileid, 1, 'uint16');
for n = 1 : EntryCount
%%%%%%%% STILL NEED TO IMPLEMENT ORIENTATION
        TagID = fread( fileid, 1, 'uint16');
	switch TagID
	 case 255, SubFileType = ReadEntry( fileid)
	 case 256, ImageWidth = ReadEntry( fileid);
	 case 257, ImageLength = ReadEntry( fileid);
	 case 258, BitsPerSample = ReadEntry( fileid);
	 case 259, Compression = ReadEntry( fileid);
	           if( Compression ~= 1), error('Cannot handle compressed images'); end;
	 case 262, PhotometricInterpretation = ReadEntry( fileid);
	 case 266, FillOrder = ReadEntry( fileid);
		   if( FillOrder ~= 1)
			if( FillOrder == 2), error('Cannot handle reverse FillOrder');
                        else error('Invalid value for FillOrder');
                        end;
 		   end
	 case 269, DocumentName = ReadEntry( fileid);
	 case 273, StripOffsets = ReadEntry( fileid);
	 case 277, SamplesPerPixel = ReadEntry( fileid);
	 case 278, RowsPerStrip = ReadEntry( fileid);
	 case 279, StripByteCounts = ReadEntry( fileid);
	 case 284, PlanarConfiguration = ReadEntry( fileid);
	           if( PlanarConfiguration < 1 | PlanarConfiguration > 2)
			 error('unhandled value for PlanarConfiguration');
		   end;
	 case 320, if( length(BitsPerSample) > 1), error('Palette images cannot have multiple BitsPerSample values'); end;
		   ColorMap = ReadEntry( fileid);
		   if( length(ColorMap) ~= 3*2^BitsPerSample), error('ColorMap has the wrong length'); end;
	 %%% IGNORE ALL OTHER TAGS
	 otherwise, fseek( fileid, 10, 0); % Skip 10 bytes
        end
end
%%% Read the offset of the next IFD
NextIFD_offset = fread( fileid, 1, 'uint32');

%%% The following tags must always be present
if(	isempty( ImageWidth)  | ...
	isempty( ImageLength)  | ...
	isempty( BitsPerSample)  | ...
	isempty( Compression)  | ...
	isempty( PhotometricInterpretation)  | ...
	isempty( StripOffsets)  | ...
	isempty( StripByteCounts)  ), error('File seems corrupted');
end

%%% DECIDE ON THE IMAGE TYPE
switch PhotometricInterpretation
 case 0, Photometric = 'Binary';
 case 1, Photometric = 'Gray';
 case 2, Photometric = 'RGB';
 case 3, if( isempty(ColorMap)), error('Colormap is missing from the file'); end;
         Photometric = 'Palette';
	 %%% LOAD COLORMAP
	 %ColorMapLength = 2^BitsPerSample;
	 %fseek( fileid, ColorMapOffset, -1);
	 %ColorMap = fread( fileid, [ColorMapLength 3], 'uint16');	 
 otherwise, error('Unhandled type of photometric interpretation');
end

%%% Construct the IFD structure to return
IFD = struct( 'ImageWidth', ImageWidth, ...
	      'ImageLength', ImageLength, ...
	      'BitsPerSample', BitsPerSample, ...
	      'SamplesPerPixel', SamplesPerPixel, ...
	      'ImageType', Photometric, ...
	      'PlanarConfig', PlanarConfiguration, ...
	      'DocumentName', DocumentName, ...
	      'StripOffsets', StripOffsets, ...
	      'StripByteCounts', StripByteCounts, ...
	      'NextIFDOffset', NextIFD_offset ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Data = ReadEntry( f) % Reads data from a single IFD entry

DataType = fread( f, 1, 'uint16');
DataCount = fread( f, 1, 'uint32');
if( DataCount < 1) error('File is corrupted (DataCount < 1)'); end;
if( DataCount == 1)
	switch DataType
	 case 2, Data = 0; %%% IGNORE 1 BYTE ASCII DATA
	 case 3, Data = fread( f, 1, 'uint16'); fseek( f, 2, 0); %SHORT
	 case 4, Data = fread( f, 1, 'uint32'); %LONG
	 otherwise, error( ['unhandled datatype: ' sprintf('%d',DataType) ]);
	end
else
	DataOffset = fread( f, 1, 'uint32');
        StoreOffset = ftell( f);
        fseek( f, DataOffset, -1);
        switch DataType
	 case 2, Data = char(fread( f, DataCount, 'uint8'))'; % ASCII
	 case 3, Data = fread( f, DataCount, 'uint16')'; % SHORT
	 case 4, Data = fread( f, DataCount, 'uint32')'; % LONG
	 otherwise, error( ['unhandled datatype: ' sprintf('%d',DataType) ]);
        end
        fseek( f, StoreOffset, -1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChanNo = find_channel_no( DocumentName)
  % Find the channel number in the Document Name

    if( isempty(DocumentName))
      ChanNo = -1;
      return;
    end;
    Keyword = 'CHANNEL';
    ChannelIdx = findstr(Keyword,upper(DocumentName));
    StrLen = length(DocumentName);
    ChanNo = getnumber( DocumentName(ChannelIdx+length(Keyword):StrLen));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function number = getnumber( string) % Gets the first number that occurs in the string

str = double(string);
L = length(str);
% Find start of number
for i = 1 : L
  if( str(i) >= '0' & str(i) <= '9') break; end;
end
% If no number was found
if( i == L)
  number = -1;
  return;
end
% Find the end of the number
start_idx = i;
for i = start_idx+1 : L
  if( ~(str(i) >= '0' & str(i) <= '9')) break; end;
end
end_idx = i-1;
% Get the actual number
numberstr = string(start_idx:end_idx);
number = sscanf(numberstr,'%i');
%-----------------------THE END ---------------------------%
