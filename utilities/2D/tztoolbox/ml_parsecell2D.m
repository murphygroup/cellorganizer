function cellcode = ...
    ml_parsecell2D( cellcode, cellbody, nucbody, da, imgsize, selcodes, param )
%ML_PARSECELL Code cell morphology.
%   CELLCODE =
%       ML_PARSECELL(CELLCODE,CELLBODY,NUCBODY,DA,IMGSIZE,SELCODES,ISSHOW)
%   returns a structure describing a cell. The input argument CELLCODE is
%   the description from previous calculation. An existing field will not
%   be updated. CELLBODY is an array of points (n1x2(3) matrix). NUCBODY is
%   also an array of points (n2x2(3) matrix). DA is the step size of coding
%   angles. IMGSIZE is the size of the original image. If ISSHOW is true,
%   the image will be shown. SELCODES is a cell array of strings, which are
%   feature names for calculation:
%       {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
%         'nuccontour','nucellhitpts','nucdist','nucelldist',...
%         'nucecc','cellarea','cellcenter','cellmangle',...
%         'cellcontour','cellhitpts','celldist','cellecc'}
%   If SELCODES is empty, all features will be submitted for calculation.
%
%   See also ML_COMBCCFEATS

% Ting Zhao
%
% Copyright (C) 2007-2013 Murphy Lab
% Carnegie Mellon University
%
% May 13, 2013 I. Cao-Berg
%    Included param structure in method so that we can pass debug and
%    display flags
%    Wrap imshow around debugging and display flags, method will save plots
%    to temporary results
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

%grj 14/5/2015
%added code to pass "display" and "debug" to updates2 block
%
%grj 17/5/2013
%fixed bug where cell and nuclear distances were only being calculated on
%the first angle
%
%icaoberg 13/1/2015
%added check for field in parameter structure

%icaoberg 13/5/2013
if nargin < 7
    error('Exactly 7 arguments are required')
end

%s2 is a structure to record which variables are ready
s2.cellbody=cellbody;
s2.nucbody=nucbody;
s2.da=da;
s2.imgsize=imgsize;

if isempty(cellcode)
    clear cellcode;
else
    %copy fields in cellcode to s2
    f=fieldnames(cellcode);
    for i=1:length(f)
        s2=setfield(s2,f{i},getfield(cellcode,f{i}));
    end
end

if isempty(selcodes)
    selcodes = {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
        'nuccontour','nucellhitpts','nucdist','nucelldist',...
        'nucecc','cellarea','cellcenter','cellmangle',...
        'cellcontour','cellhitpts','celldist','cellecc'};
end

%Calculate all required fields
s2=updates2(s2,selcodes, param );

%copy fields in s2 to cellcode
cellcode.da=s2.da;
for i=1:length(selcodes)
    cellcode=setfield(cellcode,selcodes{i},getfield(s2,selcodes{i}));
end

if param.display
    s2=updates2(s2,{'nucedge','celledge','nuccenter','len'}, param );
    img=double(s2.nucedge)+double(s2.celledge);
    for a=0:45:359
        img=ml_setimglnpixel2(img,s2.nuccenter(1:2),a,s2.len);
    end
    
    try
        imshow(img,[]);
    catch
        disp( 'Unable to display image' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s2=updates2(s2,f, param )

for i=1:length(f)
    if ~isfield(s2,f{i})
        switch f{i}
            case 'len'
                s2.len=sqrt(sum(s2.imgsize.^2));
            case 'cellarea'
                s2.cellarea=size(s2.cellbody,1);
            case 'nucarea'
                s2.nucarea=size(s2.nucbody,1);
            case 'cellecc'
                s2=updates2(s2,{'cellcenter'}, param );
                cellbody(:,1:2)=ml_addrow(s2.cellbody(:,1:2),-s2.cellcenter);
                covmat=cov(cellbody(:,1),cellbody(:,2),1);
                mu20=covmat(1,1);
                mu02=covmat(2,2);
                mu11=covmat(1,2);
                mu00=size(cellbody,1);
                medresult=sqrt((mu20-mu02)^2+4*mu11^2);
                semimajor=sqrt(2*(mu20+mu02+medresult)/mu00);
                semiminor=sqrt(2*(mu20+mu02-medresult)/mu00);
                s2.cellecc=sqrt(semimajor^2-semiminor^2)/semimajor;
            case 'nucecc'
                s2=updates2(s2,{'nuccenter'}, param );
                nucbody(:,1:2)=ml_addrow(s2.nucbody(:,1:2),-s2.nuccenter);
                covmat=cov(nucbody(:,1),nucbody(:,2),1);
                mu20=covmat(1,1);
                mu02=covmat(2,2);
                mu11=covmat(1,2);
                mu00=size(nucbody,1);
                medresult=sqrt((mu20-mu02)^2+4*mu11^2);
                semimajor=sqrt(2*(mu20+mu02+medresult)/mu00);
                semiminor=sqrt(2*(mu20+mu02-medresult)/mu00);
                s2.nucecc=sqrt(semimajor^2-semiminor^2)/semimajor;
            case 'cellimg'
                s2.cellimg=ml_obj2img2D(s2.cellbody,s2.imgsize,{'2d','bn'});
            case 'nucimg'
                s2.nucimg=ml_obj2img2D(s2.nucbody,s2.imgsize,{'2d','bn'});
            case 'celledge'
                disp('Computing cell edge');
                s2=updates2(s2,{'cellimg'}, param );
                s2.celledge=bwperim(s2.cellimg,8);
                s2.celledge=imdilate(s2.celledge, strel('disk',1));
            case 'cellcenter'
                s2.cellcenter=round(mean(s2.cellbody,1));
                %t+{27-Jan-2006}
                s2.cellcenter = s2.cellcenter(1:2);
                %t++
            case 'nuccenter'
                s2.nuccenter=round(mean(s2.nucbody,1));
                %t+{27-Jan-2006}
                s2.nuccenter = s2.nuccenter(1:2);
                %t++
            case 'nucedge'
                disp('Computing nuclear edge');
                s2=updates2(s2,{'nucimg'}, param );
                s2.nucedge=imdilate(s2.nucimg, strel('disk',1));
                s2.nucedge=bwperim(s2.nucimg,8);
            case 'cellcontour'
                disp('Computing cell contour');
                s2=updates2(s2,{'celledge'}, param );
                s2.cellcontour=ml_tracecontour(ml_mainobjcontour(s2.celledge));
            case 'nuccontour'
                disp('Computing nuclear contour');
                s2=updates2(s2,{'nucedge'}, param );
                s2.nuccontour=ml_tracecontour(ml_mainobjcontour(s2.nucedge));
            case 'cellmangle'
                s2=updates2(s2,{'cellimg'}, param );
                s2.cellmangle=ml_bwmajorangle2D(s2.cellimg)*180/pi;
            case 'nucmangle'
                s2=updates2(s2,{'nucimg'}, param );
                s2.nucmangle=ml_bwmajorangle2D(s2.nucimg)*180/pi;
            case {'celldist','cellhitpts'}
                s2=updates2(s2,{'cellcenter','len','celledge'}, param );
                s2.celldist=[];
                s2.cellhitpts=[];
                da=s2.da;
                
                if ~isempty(s2.celledge)
                    %icaoberg 20/5/2013
                    %these are neccesary when debug and display flags are
                    %set to true
                    if param.debug && param.display
                        figure
                        img = flipdim( imrotate( double(s2.celledge), 270 ), 2 );
                        imshow(img);
                        hold on
                        center = s2.cellcenter(1:2);
                        plot( center(1,1), center(1,2), 'xy' );
                    end
                    
                    for a=0:da:360-da
                        pts = ml_getlinept2(s2.cellcenter(1:2),a,s2.len);
                        ps = ml_imgptspixel(s2.celledge,pts);
                        intc2=find(ps>0);
                        
                        if ~isempty(intc2)
                            s2.cellhitpts=[s2.cellhitpts;pts(intc2(1),:)];
                            s2.celldist=[s2.celldist, ...
                                sqrt(sum((s2.cellcenter-s2.cellhitpts(end,:)).^2))];
                            
                            %icaoberg 5/20/2013
                            if param.debug && param.display
                                hold on
                                mark = pts( intc2(end), : );
                                plot( mark(1,1), mark(1,2), 'ro' )
                                if a == 0
                                    title( param.cell_image_path );
                                    ylim = get(gca,'YLim');
                                    xlim = get(gca,'XLim');
                                    angle_text = text( xlim(1), ylim(2), ['Angle:' num2str(a)], ...
                                        'VerticalAlignment', 'bottom', ...
                                        'HorizontalAlignment','left', ...
                                        'Color', 'white' );
                                else
                                    set( angle_text, ...
                                        'String', ['Angle:' num2str(a)] );
                                end
                                
                                if param.save_helper_figures
                                    F(length(F)+1) = getframe( gcf );
                                end
                                pause( param.time_to_pause );
                            end
                        else
                            %icaoberg 5/13/2013
                            if param.debug && param.display
                                figure
                                imshow(h(double(s2.celledge), ...
                                    s2.cellcenter(1:2),a,s2.len),[]);
                                title(['Not finding a hit in celldist:(image index, angle, angle step, length)=(' ...
                                    param.cell_image_filename ',' num2str(a) ',' ...
                                    num2str(da) ',' num2str(s2.len) ')'])
                                folder = [pwd filesep 'errata' ];
                                if ~exist( folder )
                                    mkdir( folder )
                                end
                                
                                filename = [folder filesep 'celldist_cell' num2str(image_index) '.jpg' ];
                                saveas( gcf, filename, 'jpg' );
                                close
                            end
                        end
                    end
                end
                close all
                
            case {'nucdist','nucelldist','nuchitpts','nucellhitpts'}
                s2=updates2(s2,{'nuccenter','len','nucedge','celledge'}, param);
                s2.nucdist=[];
                s2.nucelldist=[];
                s2.nucellhitpts=[];
                s2.nuchitpts=[];
                da=s2.da;
                
                if param.debug && param.display
                    figure
                    img = flipdim( imrotate( double(s2.nucedge), 270 ), 2 );
                    imshow(img);
                    hold on
                    center = s2.nuccenter(1:2);
                    plot( center(1,1), center(1,2), 'xy' );
                end
                
                if ~isempty(s2.nucbody)
                    for a=0:da:360-da
                        pts=ml_getlinept2(s2.nuccenter(1:2),a,s2.len);
                        ps=ml_imgptspixel(s2.nucedge,pts);
                        intc1=find(ps>0);
                        
                        if isempty(s2.cellbody)
                            intc2=1;    %just not empty
                        else
                            ps=ml_imgptspixel(s2.celledge,pts);
                            intc2=find(ps>0);
                        end
                        
                        if ~isempty(intc1) && ~isempty(intc2)
                            s2.nuchitpts=[s2.nuchitpts;pts(intc1(1),:)];
                            
                            %grj 17/5/2013
                            %set s2.nuchitpts(1,:) to s2.nuchitpts(end,:)
                            s2.nucdist=[s2.nucdist, ...
                                sqrt(sum((s2.nuccenter-s2.nuchitpts(end,:)).^2))];
                            
                            if ~isempty(s2.cellbody)
                                s2.nucellhitpts=[s2.nucellhitpts;pts(intc2(1),:)];
                                
                                %grj 17/5/2013
                                %set s2.nucellhitpts(1,:) to s2.nucellhitpts(end,:)
                                s2.nucelldist=[s2.nucelldist, ...
                                    sqrt(sum((s2.nuccenter-s2.nucellhitpts(end,:)).^2))];
                            end
                            
                            if param.debug && param.display
                                hold on
                                mark = pts( intc1(end), : );
                                plot( mark(1,1), mark(1,2), 'ro' )
                                if a == 0
                                    title( param.dna_image_path );
                                    ylim = get(gca,'YLim');
                                    xlim = get(gca,'XLim');
                                    angle_text = text( xlim(1), ylim(2), ['Angle:' num2str(a)], ...
                                        'VerticalAlignment', 'bottom', ...
                                        'HorizontalAlignment','left', ...
                                        'Color', 'white' );
                                else
                                    set( angle_text, ...
                                        'String', ['Angle:' num2str(a)] );
                                end
                                
                                pause( param.time_to_pause );
                            end
                        else
                            if isempty(s2.celledge)
                                showimg=double(s2.nucedge);
                            else
                                showimg=double(s2.nucedge)+double(s2.celledge);
                            end
                        end
                    end
                end
                close all
        end
    end
end
