% packableObject.m
%--------------------------------------------------------------------------
% This file is an example of a class which is compatible with CompressLib
% because it implements a method called "toByteArray", which outputs a 
% byte-array representation of the object which can be reconstructed
% using its constructor of the form:
%           obj = packableObject(byteArray,'PackedBytes');
%
% Objects of this class are used to test the Object-compressing 
% functionality of CompressLib in CompressLib.test.
%
% Original Author: Jesse Hopkins
%
% $Author: jhopkin $
% $Date: 2009/10/29 00:03:07 $
% $Name:  $
% $Revision: 1.2 $
%--------------------------------------------------------------------------

% $Log: packableObject.m,v $
% Revision 1.2  2009/10/29 00:03:07  jhopkin
% updated to class-version of compression methods
%

classdef packableObject
	properties
		a;
		b;
		c;
	end
	
	
	methods
		%constructor
		function obj = packableObject(varargin)
			%check inputs
			if nargin == 2 && strcmpi(varargin{2},'PackedBytes')
				%input is a packed byte-array
				byteArray = varargin{1};
				
				if ~strcmpi(class(byteArray),'uint8') || ndims(byteArray) > 2 || min(size(byteArray) ~= 1)
					error('Input must be a 1-D array of uint8');
				end
				
				%unpack the byte array to the object properties in the same
				%order that they were packed in method "toByteArray"
				
				curOffset = 1;
				%read the length of the first property
				dataLength = double(typecast(byteArray(curOffset:curOffset + 3),'uint32'));
				curOffset = curOffset + 4;
	
				%unpack "dataLength" bytes to the first property
				obj.a = CompressLib.unpackBytes(byteArray(curOffset:curOffset+dataLength-1));
				curOffset = curOffset + dataLength;
				
				%read the length of the second property
				dataLength = double(typecast(byteArray(curOffset:curOffset + 3),'uint32'));
				curOffset = curOffset + 4;
	
				%unpack "dataLength" bytes to the second property
				obj.b = CompressLib.unpackBytes(byteArray(curOffset:curOffset+dataLength-1));
				curOffset = curOffset + dataLength;
				
				%read the length of the third property
				dataLength = double(typecast(byteArray(curOffset:curOffset + 3),'uint32'));
				curOffset = curOffset + 4;
	
				%unpack "dataLength" bytes to the third property
				obj.c = CompressLib.unpackBytes(byteArray(curOffset:curOffset+dataLength-1));		
			elseif nargin == 3
				obj.a = varargin{1};
				obj.b = varargin{2};
				obj.c = varargin{3};
			else
				%default values
				obj.a = [1 2 3];
				obj.b = {'string 1','string 2','string 3'};
				obj.c = 5;
			end
		end
		
		function byteArray = toByteArray(obj)
			a_packed = CompressLib.packBytes(obj.a);
			b_packed = CompressLib.packBytes(obj.b);
			c_packed = CompressLib.packBytes(obj.c);
	
			if length(a_packed) > 2^32 || length(b_packed) > 2^32 || length(c_packed) > 2^32 
				error('could not pack table to byte array because length of byte array exceeded 2^32');
			end
			
			% pack properties into byte array as follows
			% [PropertyDataByteLength  PropertyData    PropertyDataByteLength  PropertyData ... ]
			%   uint32                   uint8[]            uint32                uint8[]
			
			byteArray = [typecast(uint32(length(a_packed)),'uint8') a_packed(:)' ...
				         typecast(uint32(length(b_packed)),'uint8') b_packed(:)' ...
				         typecast(uint32(length(c_packed)),'uint8') c_packed(:)'];
			
		end
		
	end	
end