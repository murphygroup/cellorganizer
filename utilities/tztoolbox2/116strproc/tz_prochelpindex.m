function tz_prochelpindex(filepath)
%TZ_PROCHELPINDEX processes help file generated from m2html
%   TZ_PROCHELPINDEX(FILEPATH) take the full path FILEPATH of the file and
%   modify the file for documentation.
%   
%   See also M2HTML

%   05-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

global glFundir; %Toolbox directory 

if ~exist(filepath,'file')
    error(['Could not find the file ' filepath]);
end

[fid,message] = fopen(filepath,'r');

if fid == -1 %failed to open the file
    erro(message);
end

newContents = '';
toolboxTitleMark = ['>' glFundir '/' '116'];

%Read the file line by line.
while ~feof(fid)
    contentLine = fgetl(fid);
    
    contentLine = strrep(contentLine,'Matlab Index','Toolbox Index');
   
    toolboxTitleIndices = strfind(contentLine,toolboxTitleMark)+1;
    while ~isempty(toolboxTitleIndices)
        %replace all directory names into toolbox title    
        toolboxTitleEndPositions = ...
            find( contentLine( toolboxTitleIndices(1):end ) == '<') ... 
            + (toolboxTitleIndices(1)-1) -1;
        toolboxTitle = contentLine( toolboxTitleIndices(1) : ... 
            toolboxTitleEndPositions(1) );
        toolboxTitleAbbreviation = toolboxTitle ...
            ( length(toolboxTitleMark) : end );
        
        newToolboxTitle = tz_expandabbrv(toolboxTitleAbbreviation);
        if strcmp(newToolboxTitle,'unknown word')
            newToolboxTitle = toolboxTitleAbbreviation;
        end
        
        contentLine = [ contentLine( 1 : toolboxTitleIndices(1)-1 ) ...
                newToolboxTitle ...
                contentLine( toolboxTitleEndPositions(1)+1 : end ) ];
        
        toolboxTitleIndices = strfind(contentLine,toolboxTitleMark)+1;
    end
    
    newContents = tz_addstrline(newContents,contentLine);
    
end

fclose(fid);

[fid,message] = fopen(filepath,'w');
if fid == -1 %failed to open the file
    erro(message);
end

fwrite(fid,newContents);
fclose(fid);
