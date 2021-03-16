function tz_confid_writer(outputfile, both, slf_names, id)
%ML_CONFID_WRITER   Feature confidence display
%	ML_CONFID_WRITER(OUTPUTFILE, BOTH, SLF_NAMES, THE_NOW) returns the list
%	of features to the file OUTPUTFILE given the array of type BOTH
%       and SLF_NAMES in HTML format.
%       
%
%       SImEC session identified by optional serial date number THE_NOW.
%	The features are ranked by the confidence level at which a feature
%	is distinguishable.		
%

% ml_confid_writer.m written by Edward Roques

MaxArg = 4;

if (nargin < MaxArg)
     the_now = -1;
end

num = size(both,2);

fid = fopen(outputfile, 'w');

fprintf(fid, '<HTML><HEAD><TITLE>Rank features by confidence value</TITLE></HEAD>\n');
fprintf(fid, '<BODY aLink=#0077ff bgColor=#ffffff link=#118811 text=#000000 vLink=#ff0000>\n');

fprintf(fid, '<CENTER><B><U>Rank of features by confidence value');
if (id ~= -1)
     st_time = datestr(id/10000);
     fprintf(fid, strcat(' for SImEC session started on&nbsp;', st_time));
end
fprintf(fid, '.</U></B></CENTER>\n<BR>\n\n');
fprintf(fid, '<P>Univariate t-tests were performed for each feature of the two sets of images. The highest level of significance (lowest alpha value) at which the mean feature value was found to be the same for both sets was recorded. Values of 0 indicate that the feature was never found to be the same for both images sets at any level of confidence. The resulting table therefore ranks the features according to how well they distinguish the two sets of images, features at the top of the table being better able to distinguish the sets than those at the bottom.</P>\n'); 
fprintf(fid, '<P>For further information about a feature please click the name.</P>');
fprintf(fid, '<HR>');
fprintf(fid, '<TABLE align=center cellpadding=2><TBODY>\n');
fprintf(fid, '<TR>\n');
fprintf(fid,'%33s%16s\n', '<TD><B><U>Feature name</B></U></TD>', '<TD><B><U>Confidence Level</B></U></TD>');
fprintf(fid, '</TR>\n');
for i = 1:num
fprintf(fid, '<TR>\n');
if (strncmp('Z', slf_names{both(1,i)}, 1) == 1)
     tempadd = slf_names{both(1,i)};
     tempadd(findstr(tempadd, ',')) = '_';
     fprintf(fid,'%s%s%s%s%s%12.10f%s\n', '<TD><A href="http://murphylab.web.cmu.edu/services/SLF/Z_ims/', tempadd, '.html">', slf_names{both(1,i)}, '</A></TD><TD>', both(2,i), '</TD>');
else
fprintf(fid,'%s%s%s%s%s%12.10f%s\n', '<TD><A href="http://murphylab.web.cmu.edu/services/SLF/features.html#', slf_names{both(1,i)}, '">', slf_names{both(1,i)}, '</A></TD><TD>', both(2,i), '</TD>');
end
fprintf(fid, '</TR>\n');
end

fprintf(fid, '</TBODY></TABLE>');
fprintf(fid, '<BR><BR>');
fprintf(fid, '<P align=center>Thank you for using <B>SImEC</B>.</P><P align = center><B>SImEC</B> is a service of the Murphy Lab, Carnegie Mellon University, and is available at&nbsp<A href="http://murphylab.web.cmu.edu/services/SImEC/">http://murphylab.web.cmu.edu/services/SImEC</A>.</P>');
fprintf(fid, '</BODY></HTML>\n');
status = fclose(fid);
