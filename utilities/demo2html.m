function demo2html(directory,output)
% DEMO2HTML Returns projection images and html file
% @param directory, the folder that contains the images
% @param output, the directory that save the projected images and html files

% Author: Yue Yu and Ivan Cao-Berg
% Created: Summer 2012
%
% Copyright (C) 2012 Murphy Lab
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

check_out = exist(output);
if check_out ~= 7 % the output directory is not exist
    ifcreate = input('The output directory does not exist, do you want to create it? Y/N: ','s');
    if strcmp(ifcreate,'Y')
        mkdir(output)
    else
        disp ('No directory is created')
        return
    end
end
main_html = strcat('<!DOCTYPE html><html><body><h4>',strcat(directory,'</h4><ol>'));
subfolders = dir(directory);
subfolders = subfolders(3:end);
N = length(subfolders);% the number of cells
for i = 1 : N
    tmp_folder_name = subfolders(i).name;
    tmp_html = strcat('<!DOCTYPE html><html><body><h4>',strcat(tmp_folder_name,'</h4><ol>'));
    
    proj_img_folder = strcat(output,strcat('/',tmp_folder_name));
    mkdir(proj_img_folder)% make a dir for the project image
    
    index_img_folder = strcat(directory,strcat('/',tmp_folder_name));
    image_list = dir(index_img_folder);
    image_list = image_list(3:end);
    image_num = length(image_list); % the number of index images
    for j = 1 : image_num
        img = tif2img(strcat(index_img_folder,strcat('/',image_list(j).name))); 
        % the current image in current folder
        out_img = img2projection(img);
        tmp_img_name = regexp(image_list(j).name,'\.','split');
        tmp_img_name = tmp_img_name{1};
        tmp_img_name = strcat(tmp_img_name,'.jpg');
        imwrite(uint8(ml_bcimg(out_img,[],[])),strcat(proj_img_folder,strcat('/',tmp_img_name)),'jpg');
        %add the image to html list
        tmp_html = strcat(tmp_html,'<li><a');
        tmp_html = strcat(tmp_html,' href="');
        tmp_html = strcat(tmp_html,tmp_img_name);
        tmp_html = strcat(tmp_html,'">');
        tmp_html = strcat(tmp_html,tmp_img_name);
        tmp_html = strcat(tmp_html,'</a></li>');
    end
    tmp_html = strcat(tmp_html,'</ol></body></html>');
    fid = fopen(strcat(proj_img_folder,'/index.html'),'w');
    fprintf(fid,'%s',tmp_html);
    fclose(fid);
    
    %make the main html
    main_html = strcat(main_html,'<li><a');
    main_html = strcat(main_html,' href="');
    main_html = strcat(main_html,strcat('./',tmp_folder_name));
    main_html = strcat(main_html,'/index.html');
    main_html = strcat(main_html,'">');
    main_html = strcat(main_html,tmp_folder_name);
    main_html = strcat(main_html,'</a></li>');
end
main_html = strcat(main_html,'</ol></body></html>');
fid = fopen(strcat(output,'/index.html'),'w');
fprintf(fid,'%s',main_html);
fclose(fid);
