% 2011-11-15 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image.
% 2011-12-07 tebuck: Adding naming_convention.blank_channels to
% blank any channels for feature computation, not just as with
% channel_as_protein.

%if isempty(erstruct.channel)
if isempty(erstruct.channel) && ~(er_channel_blank && optimize)
    erstruct.channel = imread(erstruct.channel_path);
    %if strcmpi(channel_as_protein, 'er') || sum(strcmpi('er', naming_convention.blank_channels)) > 0
    %if strcmpi(channel_as_protein, 'er') || er_channel_blank
    if er_channel_blank
      %protstruct.channel = erstruct.channel; 
      erstruct.channel = erstruct.channel * 0; 
    end
end


