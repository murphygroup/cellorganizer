% 2011-11-15 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image.
% 2011-12-07 tebuck: Adding naming_convention.blank_channels to
% blank any channels for feature computation, not just as with
% channel_as_protein.

%if isempty(nucstruct.channel)
if isempty(nucstruct.channel) && ~(nuclear_channel_blank && optimize)
    nucstruct.channel = imread(nucstruct.channel_path);
    %if strcmpi(channel_as_protein, 'nucleus') || sum(strcmpi('nucleus', naming_convention.blank_channels)) > 0
    %if strcmpi(channel_as_protein, 'nucleus') || nuclear_channel_blank
    if nuclear_channel_blank
      %protstruct.channel = nucstruct.channel; 
      nucstruct.channel = nucstruct.channel * 0; 
    end
end

