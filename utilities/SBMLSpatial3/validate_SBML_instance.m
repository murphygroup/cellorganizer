function [sbml_valid, sbml_problems] = validate_SBML_instance( filename, save_validation_results )

% Author: Ivan E. Cao-Berg, Taraz Buck
%
% Copyright (C) 2016-2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

if ~exist(filename, 'file')
    error('file ''%s'' does not exist', filename);
end

if nargin < 2
    save_validation_results = false;
end

sbml_problem_blank = struct('category', '', 'code', '', 'severity', '', 'location_line', [], 'location_column', [], 'message_lang', '', 'message_text', '', 'package_code', '', 'shortmessage_lang', '', 'shortmessage_text', '', 'detail_category', '', 'detail_severity', '');

sbml_valid = false;
sbml_problems = repmat(sbml_problem_blank, 0, 1);

filename_with_content_type = filename;
if strcmpi(filename(end-3:end), '.xml')
    filename_with_content_type = [filename_with_content_type, ';type=application/xml'];
elseif strcmpi(filename(end-7:end), '.xml.zip')
    error('Not working yet, .xml.zip validates when .xml fails');
    % filename_with_content_type = [filename_with_content_type, ';type=application/xml+zip'];
    filename_with_content_type = [filename_with_content_type, ';type=application/zip'];
else
    error('Unknown file type');
end

command_start_time = tic();

% command = ['curl -F "file=@' filename_with_content_type '" -F output=xml -F offcheck=u' ...
    % ' --retry 2 --retry-max-time 30 http://sbml.org/validator/'];

status = 1;
attempt_count = 0;
attempt_max = 3;
attempt_timeout = 20;
while status ~= 0 && attempt_count < attempt_max
    attempt_count = attempt_count + 1;
    command = sprintf(['curl -F "file=@%s" -F output=xml -F offcheck=u' ...
        ' --max-time %i http://sbml.org/validator/'], ...
        filename_with_content_type, attempt_timeout);
    [status,result] = system(command);
end
if status ~= 0
    error('status %i after %i attempts, result ''%s''', status, attempt_count, result);
end

result = strtrim(result);

if save_validation_results
    fid = fopen([filename, '.validationresults.xml'], 'w');
    fwrite(fid, result);
    fclose(fid);
end

% disp( result );
command_elapsed_time = toc(command_start_time);
fprintf('\n%s command took %.3f s\n\n', mfilename, command_elapsed_time);

docNode = xmlread(org.xml.sax.InputSource(java.io.StringReader(result)));
docRootNode = docNode.getDocumentElement;
error_nodes = docNode.getElementsByTagName('error');
problem_nodes = docNode.getElementsByTagName('problem');
sbml_valid = error_nodes.getLength() == 0;
sbml_problems = repmat(sbml_problem_blank, 0, 1);
for i = 1:problem_nodes.getLength()
    problem_node = problem_nodes.item(i-1);
    problem = sbml_problem_blank;
    try
        problem.category = char(problem_node.getAttribute('category'));
    end
    try
        problem.code = char(problem_node.getAttribute('code'));
    end
    try
        problem.severity = char(problem_node.getAttribute('severity'));
    end
    try
        problem.location_line = char(problem_node.getElementsByTagName('location').item(0).getAttribute('line'));
    end
    try
        problem.location_column = char(problem_node.getElementsByTagName('location').item(0).getAttribute('column'));
    end
    try
        problem.message_lang = char(problem_node.getElementsByTagName('message').item(0).getAttribute('lang'));
    end
    try
        problem.message_text = char(problem_node.getElementsByTagName('message').item(0).getTextContent());
    end
    try
        problem.package_code = char(problem_node.getElementsByTagName('package').item(0).getAttribute('code'));
    end
    try
        problem.shortmessage_lang = char(problem_node.getElementsByTagName('shortmessage').item(0).getAttribute('lang'));
    end
    try
        problem.shortmessage_text = char(problem_node.getElementsByTagName('shortmessage').item(0).getTextContent());
    end
    try
        problem.detail_category = char(problem_node.getElementsByTagName('detail').item(0).getAttribute('category'));
    end
    try
        problem.detail_severity = char(problem_node.getElementsByTagName('detail').item(0).getAttribute('severity'));
    end
    
    sbml_problems(end+1) = problem;
end
warning('Does not process warnings');

