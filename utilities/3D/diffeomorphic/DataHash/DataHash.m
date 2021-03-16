function Hash = DataHash(Data, Opt)
% DATAHASH - Create a checksum for Matlab array of any type
% Hash = DataHash(Data, Opt)
% INPUT:
%   Data: Array of these built-in types:
%           (U)INT8/16/32/64, SINGLE, DOUBLE, (real or complex)
%           CHAR, LOGICAL, CELL, STRUCT (scalar or array, nested),
%           function_handle.
%   Opt:  Struct to specify the hashing algorithm and the output format.
%         Opt and all its fields are optional.
%         Opt.Method: String, known methods for Java 1.6 (Matlab 2009a):
%              'SHA-1', 'SHA-256', 'SHA-384', 'SHA-512', 'MD2', 'MD5'.
%            Known methods for Java 1.3 (Matlab 6.5):
%              'MD5', 'SHA-1'.
%            Default: 'MD5'.
%         Opt.Format: String specifying the output format:
%            'hex', 'HEX':      Lower/uppercase hexadecimal string.
%            'double', 'uint8': Numerical vector, same values.
%            'base64':          Base64 encoded string, only printable
%                               ASCII characters, 33% shorter than 'hex'.
%             Default: 'hex'.
%         Opt.isFile: LOGICAL flag, if TRUE the 1st input is the file name the
%             hash is created for. Default: FALSE.
%
% OUTPUT:
%   Hash: String or numeric vector. The number of elements depends on the
%         hashing method.
%
% EXAMPLES:
% % Default: MD5, hex:
%   DataHash([])                % 7de5637fd217d0e44e0082f4d79b3e73
% % MD5, Base64:
%   Opt.Format = 'base64';
%   Opt.Method = 'MD5';
%   DataHash(int32(1:10), Opt)  % bKdecqzUpOrL4oxzk+cfyg
% % SHA-1, Base64:
%   S.a = uint8([]);
%   S.b = {{1:10}, struct('q', uint64(415))};
%   Opt.Method = 'SHA-1';
%   DataHash(S, Opt)            % ZMe4eUAp0G9TDrvSW0/Qc0gQ9/A
%
% NOTE:
%   Currently objetcs are not handled, but they can if the user provides a
%   function to typecast the values to a UINT8 vector. Perhaps STRUCT(Obj) is
%   sufficient.
%   Function_handles are converted to a struct by the built-in function
%   FUNCTIONS which leads to dependencies to the platform and Matlab version.
%   The behaviour can be adjusted in the subfunction ConvertFuncHandles.
%
%   This function tries to use James Tursa's TYPECASTX, because it is smarter
%   and much faster than the built-in TYPECAST:
%     http://www.mathworks.com/matlabcentral/fileexchange/17476
%   For Matlab 6.5 installing TYPECASTX is obligatory to run DataHash.
%   @James: Thanks!
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, Java: 1.3.1_01, 1.6.0_04.
% Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also: TYPECAST, CAST.
% FEX:
% Michael Kleder, "Compute Hash", no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/8944
% Tim, "Serialize/Deserialize", converts structs and cells to a byte stream:
%   http://www.mathworks.com/matlabcentral/fileexchange/29457
% Jan Simon, "CalcMD5", MD5 only, faster C-mex, no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/25921

% $JRev: R-a V:001 Sum:VHCXCEMfTsyI Date:01-May-2011 21:52:09 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Published\DataHash\DataHash.m $
% History:
% 001: 01-May-2011 21:52, First version.

% Main function: ===============================================================
% Java is needed:
if ~usejava('jvm')
   error(['JSimon:', mfilename, ':NoJava'], ...
      ['*** ', mfilename, ': Java is needed.']);
end

% James Tursa's TYPECASTX creates a shared data copy instead of the deep copy as
% Matlab's TYPECAST - for a [1000x1000] DOUBLE array this is 100 times faster!
useTYPECASTX = ~isempty(which('typecastx'));

% Default options: -------------------------------------------------------------
Method    = 'MD5';
OutFormat = 'hex';
isFile    = false;

% Check number and type of inputs:
nArg = nargin;
if nArg == 2
   if isa(Opt, 'struct')
      % Specify hash algorithm:
      if isfield(Opt, 'Method')
         Method = upper(Opt.Method);
      end
      
      % Specify output format:
      if isfield(Opt, 'Format')
         OutFormat = Opt.Format;
      end
      
      % [Data] is a file name, if [Opt.isFile] is enabled:
      if isfield(Opt, 'isFile')
         isFile = any(Opt.isFile);
         if isFile
            if ischar(Data) == 0
               error(['JSimon:', mfilename, ':CannotOpen'], ...
                  ['*** ', mfilename, ': [Opt.isFile] is enabled, ', ...
                     'but 1st input is no file name.']);
            end
            
            if exist(Data, 'file') ~= 2
               error(['JSimon:', mfilename, ':FileNotFound'], ...
                  ['*** ', mfilename, ': File not found: %s.'], Data);
            end
         end
      end
      
   else  % Bad type of 2nd input:
      error(['JSimon:', mfilename, ':BadInput2'], ...
         ['*** ', mfilename, ': 2nd input [Opt] must be a struct.']);
   end
   
elseif nArg ~= 1  % Bad number of arguments:
   error(['JSimon:', mfilename, ':BadNInput'], ...
      ['*** ', mfilename, ': 1 or 2 inputs required.']);
end

% Create the engine: -----------------------------------------------------------
try
   Engine = java.security.MessageDigest.getInstance(Method);
catch
   error(['JSimon:', mfilename, ':BadInput2'], ...
      ['*** ', mfilename, ': Invalid algorithm: [%s].'], Method);
end

% Create the hash value: -------------------------------------------------------
if isFile
   % Read the file and calculate the hash:
   FID = fopen(Data, 'r');
   if FID < 0
      error(['JSimon:', mfilename, ':CannotOpen'], ...
         ['*** ', mfilename, ': Cannot open file: %s.'], Data);
   end
   Data = fread(FID, Inf, '*uint8');
   fclose(FID);
   Engine.update(Data);
   if useTYPECASTX
      Hash = typecastx(Engine.digest, 'uint8');
   else
      Hash = typecast(Engine.digest, 'uint8');
   end
   
elseif useTYPECASTX  % Use TYPECASTX:
   Engine = CoreHash_X(Data, Engine);
   Hash   = typecastx(Engine.digest, 'uint8');
   
else                 % Use the slower built-in TYPECAST:
   Engine = CoreHash(Data, Engine);
   Hash   = typecast(Engine.digest, 'uint8');
end

% Convert hash specific output format: -----------------------------------------
switch OutFormat
   case 'hex'
      Hash = sprintf('%.2x', double(Hash));
   case 'HEX'
      Hash = sprintf('%.2X', double(Hash));
   case 'double'
      Hash = double(reshape(Hash, 1, []));
   case 'uint8'
      Hash = reshape(Hash, 1, []);
   case 'base64'
      Hash = fBase64_enc(double(Hash));
   otherwise
      error(['JSimon:', mfilename, ':BadOutFormat'], ...
         '*** %s: [Opt.Format] must be: HEX, hex, uint8, double, base64.', ...
         mfilename);
end

% return;

% ******************************************************************************
function Engine = CoreHash_X(Data, Engine)
% This mothod uses the fast TYPECASTX written by James Tursa.

% Consider the type and dimensions of the array to distinguish arrays with the
% same data, but different shape: [0 x 0] and [0 x 1], [1,2] and [1;2],
% DOUBLE(0) and SINGLE([0,0]):
Engine.update([uint8(class(Data)), typecastx(size(Data), 'uint8')]);

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash_X(F, Engine);   % Catch the fieldnames
   
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash_X(Data(iS).(F{iField}), Engine);
      end
   end
   
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash_X(Data{iS}, Engine);
   end
      
elseif isnumeric(Data) || islogical(Data) || ischar(Data)
   if isempty(Data) == 0
      if isreal(Data)                % TRUE for LOGICAL and CHAR also:
         Engine.update(typecastx(Data(:), 'uint8'));
      else                           % TYPECASTX accepts complex input:
         tmp = typecastx(Data(:), 'uint8');
         Engine.update(real(tmp));
         Engine.update(imag(tmp));
      end
   end
   
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
   
else
   % If the user can provide a conversion to a UINT8 vector, considering objects
   % would be easy...
   warning(['JSimon:', mfilename, ':BadDataType'], ...
      ['Type of variable not considered: ', class(Data)]);
end

% return;

% ******************************************************************************
function Engine = CoreHash(Data, Engine)
% This methods uses the slower TYPECAST of Matlab
% See CoreHash_X for comments.

Engine.update([uint8(class(Data)), typecast(size(Data), 'uint8')]);

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash(F, Engine);     % Catch the fieldnames
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash(Data(iS).(F{iField}), Engine);
      end
   end
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash(Data{iS}, Engine);
   end
elseif isempty(Data)
elseif isnumeric(Data)
   if isreal(Data)
      Engine.update(typecast(Data(:), 'uint8'));
   else
      Engine.update(typecast(real(Data(:)), 'uint8'));
      Engine.update(typecast(imag(Data(:)), 'uint8'));
   end
elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
elseif ischar(Data)                  % TYPECAST cannot handle CHAR
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
else
   warning(['JSimon:', mfilename, ':BadDataType'], ...
      ['Type of variable not considered: ', class(Data)]);
end

% return;

% ******************************************************************************
function FuncKey = ConvertFuncHandle(FuncH)

%   The subfunction ConvertFuncHandle converts function_handles to a struct
%   using the Matlab function FUNCTIONS. The output of this function changes
%   with the Matlab version, such that DataHash(@sin) replies different hashes
%   under Matlab 6.5 and 2009a.
%   An alternative is using the function name and name of the file for
%   function_handles, but this is not unique for nested or anonymous functions.
%   If the MATLABROOT is removed from the file's path, at least the hash of
%   Matlab's toolbox functions is (usually!) not influenced by the version.
%   Finally I'm in doubt if there is a unique method to hash function handles.
%   Please adjust the subfunction ConvertFuncHandles to your needs.

% The Matlab version influences the conversion by FUNCTIONS:
% 1. The format of the struct replied FUNCTIONS is not fixed,
% 2. The full paths of toolbox function e.g. for @mean differ.
FuncKey = functions(FuncH);

% ALTERNATIVE: Use name and path. The <matlabroot> part of the toolbox functions
% is replaced such that the hash for @mean does not depend on the Matlab
% version.
% Drawbacks: Anonymous functions, nested functions...
% funcStruct = functions(FuncH);
% funcfile   = strrep(funcStruct.file, matlabroot, '<MATLAB>');
% FuncKey    = uint8([funcStruct.function, ' ', funcfile]);

% Finally I'm afraid there is no unique method to get a hash for a function
% handle. Please adjust this conversion to your needs.

% return;

% ******************************************************************************
function Out = fBase64_enc(In)
% Encode numeric vector of UINT8 values to base64 string.

Pool = [65:90, 97:122, 48:57, 43, 47];  % [0:9, a:z, A:Z, +, /]
v8   = [128; 64; 32; 16; 8; 4; 2; 1];
v6   = [32, 16, 8, 4, 2, 1];

In  = reshape(In, 1, []);
X   = rem(floor(In(ones(8, 1), :) ./ v8(:, ones(length(In), 1))), 2);
Y   = reshape([X(:); zeros(6 - rem(numel(X), 6), 1)], 6, []);
Out = char(Pool(1 + v6 * Y));

% return;
