function [LeftX,LeftY,RightX,RightY]=ExchangeExcelToOur(tmp_cell)
% Converts the point data in Excel to the left and right
LeftX=zeros(size(tmp_cell,1),1);
LeftY=zeros(size(tmp_cell,1),1);
RightX=zeros(size(tmp_cell,1),1);
RightY=zeros(size(tmp_cell,1),1);
for k=1:size(tmp_cell,1)
    distance=tmp_cell(k,1);
    Angle=tmp_cell(k,2);
    Left=tmp_cell(k,3);
    Top=tmp_cell(k,4);
    Width=tmp_cell(k,5);
    Height=tmp_cell(k,6);
    angle=Angle*pi/180;% Turn the angle to radians
    x=Left+Width/2;
    y=Top+Height/2;
    seg_offset(1)=cos(-angle)*distance*0.5;
    seg_offset(2)=sin(-angle)*distance*0.5;
    LeftX(k,1)=x-seg_offset(1);
    LeftY(k,1)=y-seg_offset(2);
    RightX(k,1)=x+seg_offset(1);
    RightY(k,1)=y+seg_offset(2);
end