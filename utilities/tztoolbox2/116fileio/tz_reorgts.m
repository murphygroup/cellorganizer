function tz_reorgts(tsdir)
%TZ_REORGTS Reorganize time series images from QED.
%   TZ_REORGTS(TSDIR) separates the time series images under TSDIR into
%   several directories, each of which contains the images at one time
%   point. This might not work for the new version of QED.
%   
%   See also

%   ??-???-???? Initial write T. Zhao
%   14-DEC-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

image_list = ml_dir([tsdir '/' '*.tif']);
for image_no = 1:length(image_list)
    image_name=image_list{image_no};
    time_points(image_no)=str2num(image_name(end-4))
end

time_num=max(time_points);

for time_point=1:time_num
    tpdir=[tsdir '-000'];
    strindex=num2str(time_point);
    tpdir(end-length(time_point)+1:end)=strindex;
    if(~exist(tpdir,'dir'))
        tz_make_dir(tpdir);
    end
%    L=length(tpdir);
    tpdirlist(time_point)={tpdir};
end
 
for image_no = 1:length(image_list)
    image_list{image_no}
    tpdirlist{time_points(image_no)}
    command=['cp ' tsdir '/' image_list{image_no} ' ' tpdirlist{time_points(image_no)}];
    unix(command);
end
