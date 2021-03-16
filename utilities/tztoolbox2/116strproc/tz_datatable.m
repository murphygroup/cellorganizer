function tablestr = tz_datatable(x,rownames,colnames,format)
%TZ_DATATABLE Generate a string of data table.
%   TABLESTR = TZ_DATATABLE(X,ROWNAMES,COLNAMES,FORMAT) returns a string
%   for displaying data X in a table. X is a numeric matrix. ROWNAMES 
%   is a [string array] of names
%   for each row. COLNAMES is a [string array] of names for each column.
%   ROWNAMES or COLNAMES is empty if there are no such names. FORMAT is 
%   a string specifying the format of the table:
%       'html' - html. X can be a cell array.
%       'latex' - latex
%       'common','common2' - for matlab displaying
%   
%   See also

%   31-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 4
    error('Exactly 4 arguments are required')
end

lineBreak = char(10);

switch format
    case 'html'
        tablestr = ['<table border=1>' lineBreak];
        
        if ~isempty(colnames)
            tablestr = [tablestr '<tr>'];
            if ~isempty(rownames)
                tablestr = [tablestr '<td></td>'];
            end
            for i=1:length(colnames)
                tablestr = [tablestr '<td>' colnames{i} '</td>'];
            end
            tablestr = [tablestr '</tr>' lineBreak];
        end
        
        for i=1:size(x,1)
            tablestr = [tablestr '<tr>'];
            if ~isempty(rownames)
                tablestr = [tablestr '<td>' rownames{i} '</td>'];
            end
            for j=1:size(x,2)
                if iscell(x)
                    cellstr = num2str(x{i,j});
                else                                   
                    if i==j                                       
                        cellstr = ['<b>' num2str(x(i,j)) '</b>'];
                    else
                        cellstr = num2str(x(i,j));
                    end
                end

                tablestr = [tablestr '<td>' cellstr '</td>'];
            end
            tablestr = [tablestr '</tr>' lineBreak];
        end
        tablestr = [tablestr '</table>' lineBreak];
    case 'latex'
        tablestr = ['\begin{table}' lineBreak];
        tablestr = [tablestr '\begin{center}' lineBreak];
        
        %\begin{tabular}{c|cccccc...}
        tablestr = [tablestr '\begin{tabular}'];
        tablestr = [tablestr '{'];
        
        if ~isempty(rownames)
            tablestr = [tablestr 'c|'];
        end
        
        for i=1:size(x,2)
           tablestr = [tablestr 'c']; 
        end
        tablestr = [tablestr '}'];
        
        tablestr = [lineBreak tablestr '\hline' lineBreak];
        
        if ~isempty(colnames)
            tablestr = ...
                [tablestr '&' tz_cell2str(colnames,'&') '\\' lineBreak];
        end
        tablestr = [tablestr '\hline' lineBreak];
        
        for i=1:size(x,1)
            if ~isempty(rownames)
                tablestr = [tablestr rownames{i}];
            end
            for j=1:size(x,2)
                tablestr = [tablestr '&' num2str(x(i,j))];
            end
            tablestr = [tablestr '\\' lineBreak];
        end
        
        tablestr = [tablestr '\hline' lineBreak];
        tablestr = [tablestr '\end{tabular}' lineBreak];
        tablestr = [tablestr '\end{center}' lineBreak];
        tablestr = [tablestr '\end{table}' lineBreak];
    case 'common'
        tablestr = ['Rows:' tz_cell2str(rownames) lineBreak];
        tablestr = [tablestr 'Columns:' tz_cell2str(colnames) lineBreak];
        tablestr = [tablestr lineBreak];
        for i=1:size(x,1)
            for j=1:size(x,2)
                tablestr = [sprintf([tablestr '%8.2e '],x(i,j))];
            end
            tablestr = [tablestr lineBreak];
        end
    case 'common2'
        tablestr = ['Rows:' tz_cell2str(rownames) lineBreak];
        tablestr = [tablestr 'Columns:' tz_cell2str(colnames) lineBreak];
        tablestr = [tablestr lineBreak];
        for i=1:size(x,1)
            tablestr = [tablestr num2str(x(i,:)) lineBreak];
        end
    otherwise
        error('Unrecognized format name.');
end
