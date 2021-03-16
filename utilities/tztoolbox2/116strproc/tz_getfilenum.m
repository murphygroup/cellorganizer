function num = tz_getfilenum(filename,pos)

%TZ_GETFILENUM Extract number label from a file name.
%   TZ_GETFILENUM(FILENAME) extracts a number of the last occurrence from
%   FILENAME.
%
%   TZ_GETFILENUM(FILENAME,POS) extracts a number from FILENAME. The
%   occurrence of the number is specified by POS. If POS=1, it is the
%   last occurrence.
%
%   Example:
%   tz_getfilenum('image0403-1-023.tif') returns 23.
%   tz_getfilenum('image0403-1-023.tif',2) returns 1.
%   tz_getfilenum('image0403-1-023.tif',3) returns 403.

%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%   15-MAR-2004 Modified TINGZ
%       - add pos
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_getfilenum','ml_getfilenum'));

if nargin < 2
   pos = 1;
end

filename = fliplr(filename);
num = [];
curpos = 1;

for i = 1:length(filename)
   if filename(i) >= '0' & filename(i) <= '9'
      num = [num filename(i)];
   else
      if ~isempty(num)
         if curpos == pos
            break
         else
            curpos = curpos+1;
            num = [];
         end
      end
   end
end

if ~isempty(num)
   num = fliplr(num);
   num = str2num(num);
else
   warning('no number found');
   num = -1;
end
