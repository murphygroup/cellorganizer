
%if isempty(protstruct.channel)
if isempty(protstruct.channel) && ~(protein_channel_blank && optimize)
    protstruct.channel = imread(protstruct.channel_path);
% $$$     % Optionally switch a channel with the protein channel:
% $$$     switch lower(channel_as_protein)
% $$$      case 'nuclear'
% $$$       protstruct.channel = imread(nucstruct.channel_path);
% $$$      case 'tubulin'
% $$$       protstruct.channel = imread(tubstruct.channel_path);
% $$$      case 'er'
% $$$       protstruct.channel = imread(erstruct.channel_path);
% $$$      case 'protein'
% $$$       protstruct.channel = imread(protstruct.channel_path);
% $$$     end
    if protein_channel_blank
      protstruct.channel = protstruct.channel * 0; 
    end
end
