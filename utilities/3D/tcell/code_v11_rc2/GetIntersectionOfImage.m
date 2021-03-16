function [PointIdex,GfpOutPut,DICOutPut,OverlayOutPut]=GetIntersectionOfImage(image1,image2,image3,FrameNum,Width)
% image1 GFP of initial image
% image2 DIC of initial image
% image3 Overlay of initial image

% PointIdex is an array of structures  containing intersection data  and fitting data of the corresponding big cell and small cell
% Width is the width of the deducted small cell
% FrameNum  is   number of frames to be processed

% segmentation small cells using the method of Watershed
loc=MarkerControlled_Watershed2(image1);
%  figure(FrameNum);clf;imshow(image2);

% detection of large cells
[centers, radii]=DrawBigCell(image2,15,35);
% preserve large cells and small cells with intersecting points  into the struct of PointIdex
PointIdex=struct();
cnt=1;
for k=1:length(loc)
    % hold on;plot(loc{k,1}(1),loc{k,1}(2),'g*');hold on;   % draw the centroid of the detected small cell
    SmallCellData=Process_SmallCell(image1,[loc{k,1}(1),loc{k,1}(2)],Width,Width);
    % Save the point about small cell and the big cell has the point of intersection
    if isfield(SmallCellData,'Zhou')
        [JiaoDian_Point,IdexBigCell]=GetNode(SmallCellData,radii,centers);
        if ~isempty(IdexBigCell)
            PointIdex(cnt).BigCellIdex=IdexBigCell;
            PointIdex(cnt).BigCellCenter=centers(IdexBigCell,:);   % Save the center of the big cell
            PointIdex(cnt).BigCellRadii=radii(IdexBigCell);        % Save the radius of the big cell
            PointIdex(cnt).SmallCenter=loc{k,1};                   % Save the centroid position of small cell
            PointIdex(cnt).SmallCellIdex=k;                        % Save the number of the small cell with the intersection
            PointIdex(cnt).SmallCellData=SmallCellData;            % Save the relevant information of the small cell fit
            PointIdex(cnt).JiaoDian_Point=JiaoDian_Point;          % Save the intersection position                      
            cnt=cnt+1;
        end
    end
end
% Show contours and intersections
% for h=1:cnt-1
%     hold on;plot(PointIdex(h).SmallCellData.Center(1),PointIdex(h).SmallCellData.Center(2),'r*');% Draw the center of the ellipse with the intersection
%     hold on;plot(centers(PointIdex(h).BigCellIdex,1),centers(PointIdex(h).BigCellIdex,2),'b*');% Draw the center of the big circle with the intersection
%     hold on;viscircles(centers(PointIdex(h).BigCellIdex,:), radii(PointIdex(h).BigCellIdex),'EdgeColor','b');% Draw a large circle with intersections
%     % Obtain small cell ellipse parameters
%     p=PointIdex(h).SmallCellData.equ{2};F=PointIdex(h).SmallCellData.equ{4};
%     xmin=PointIdex(h).SmallCellData.mimaxXY(1);xmax=PointIdex(h).SmallCellData.mimaxXY(2);
%     ymin=PointIdex(h).SmallCellData.mimaxXY(3);ymax=PointIdex(h).SmallCellData.mimaxXY(4);
%     hold on;h1=ezplot(@(x,y)F(p,[x,y]),[xmin,xmax,ymin,ymax]);set(h1,'Color','r','LineWidth',3);  % Draw an ellipse
%     hold on;plot(PointIdex(h).JiaoDian_Point(:,1),PointIdex(h).JiaoDian_Point(:,2),'k.','MarkerSize',25);
% end