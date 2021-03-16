function [ Forward_Need,Need] = No_Touch( Start,Overlay,DIC,GFP,Excel_lp, Excel_rp,OnsetTime,Need,Processdir,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 if Start > 2
       [Forward_Need]=TurnFront2(Overlay,DIC,GFP, Start,Excel_lp, Excel_rp,Processdir,param);
       Forward_Need=[Forward_Need,OnsetTime * ones(size(Forward_Need, 1), 1)];    
%         if ~isempty(Forward_Need)
%             if  Forward_Need(1,1)==0
%                 initpoint=Forward_Need(2,:);
%             else
%                 initpoint=Forward_Need(1,:);
%             end
%             figure(1)
%             imshow(DIC{initpoint(5),1});
%             hold on
%             plot(initpoint(1),initpoint(2),'r*',initpoint(3),initpoint(4),'r*');
%             JD  = Det_initpoint( initpoint,DIC,GFP,Start);
%            %判断是否需要使用这个点当做初始点
%             if  (abs(JD(1)+JD(3)-Need(1)-Need(3))+abs(JD(2)+JD(4)-Need(2)-Need(4)) )/4 > 5 || JD(5)~=Need(5)
%                 Start=length(GFP);
%                 Forward_Need = [];
%                 Need = [];
%             end
%         end
    elseif Start==2
       [Forward_Need]=TurnFront2(Overlay,DIC,GFP, Start,Excel_lp, Excel_rp,Processdir,param);
       Forward_Need=[Forward_Need,OnsetTime * ones(size(Forward_Need, 1), 1)];    
%         if ~isempty(Forward_Need)
%             if  Forward_Need(1,1)==0
%                 initpoint=Forward_Need(2,:);
%             else
%                 initpoint=Forward_Need(1,:);
%             end
%             JD  = Det_initpoint( initpoint,DIC,GFP,Start);
%           
%              if  (abs(JD(1)+JD(3)-Need(1)-Need(3))+abs(JD(2)+JD(4)-Need(2)-Need(4)) )/4 > 5 || JD(5)~=Need(5)
%                 Start=length(GFP);
%                 Forward_Need = [];
%                 Need = [];
%             end
%         end
 else
          Forward_Need=[];
    end

end

