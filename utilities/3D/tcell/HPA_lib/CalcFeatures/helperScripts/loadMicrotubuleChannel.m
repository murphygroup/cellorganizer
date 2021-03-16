% 2011-11-15 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image.
% 2011-12-07 tebuck: Adding naming_convention.blank_channels to
% blank any channels for feature computation, not just as with
% channel_as_protein.

%if isempty(tubstruct.channel)
if isempty(tubstruct.channel) && ~(tubulin_channel_blank && optimize)
    tubstruct.channel = imread(tubstruct.channel_path);
    %if strcmpi(channel_as_protein, 'tubulin') || sum(strcmpi('tubulin', naming_convention.blank_channels)) > 0
    %if strcmpi(channel_as_protein, 'tubulin') || tubulin_channel_blank
    if tubulin_channel_blank
      %protstruct.channel = tubstruct.channel; 
      tubstruct.channel = tubstruct.channel * 0; 
    end
end


