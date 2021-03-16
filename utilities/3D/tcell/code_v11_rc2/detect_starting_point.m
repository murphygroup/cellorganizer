function [SaveName] = detect_starting_point(FilePath1, SaveName )

StrDic= FilePath1;
StrOverlay=StrDic;StrGfp=StrDic;

StartFrame=5;
All_Start=[];
% Detection;All_Start variable is all detection point
GFP=GfpDeal( StrGfp);         
DIC=DICDeal( StrDic  );
Overlay=OverlayDeal( StrOverlay );
while(1)
    if StartFrame>40
        break;
    end
    image1 = GFP{StartFrame,1};
    image2 = DIC{StartFrame,1};
    image3 = Overlay{StartFrame,1};
    %The intersection point is obtained for each start image
    [PointIdex]=GetIntersectionOfImage(image1,image2,image3,StartFrame,34);
    %Each intersection of each initial frame are pushed forward to get the intersection
    for k=1:length(PointIdex)
        Point=PointIdex(k).JiaoDian_Point;
        Excel_lp=Point(1,:);Excel_rp=Point(2,:);
        Need  = TurnFront( Overlay,DIC,GFP, StartFrame,Excel_lp, Excel_rp);
        [All_Start] = GainStart( Need,All_Start );
        Need=[];
    end
    StartFrame=StartFrame+5;
end
% Save the result determined as the initial point
FormHead={'Lx','Ly','Rx','Ry','start_frame'};
All_Start=num2cell(All_Start);
Result = cell2table(All_Start,'VariableNames',FormHead);
writetable(Result,SaveName,'WriteRowNames',true);
end