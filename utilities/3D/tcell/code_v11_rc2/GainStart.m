function [All_Start] = GainStart( Need,All_Start )
m=find(Need(:,1)==0);
if ~isempty(m)
    %Only one point do not
    if length(m)+1==size(Need,1)
        All_Start=All_Start;
    else
        Temp=Need(m(1)-1,:);
        % if record this point and all the points in front of is approaching,it will do not record
        for i=1:size(All_Start,1)
            if (abs(All_Start(i,1)-Temp(1))+abs(All_Start(i,2)-Temp(2))+abs(All_Start(i,3)-Temp(3))+abs(All_Start(i,4)-Temp(4)))/4<15
                Temp=[];
                break;
            end
        end       
        All_Start=[All_Start;Temp];
    end
else
    Temp=Need(size(Need,1),:);
    for i=1:size(All_Start,1)
        if (abs(All_Start(i,1)-Temp(1))+abs(All_Start(i,2)-Temp(2))+abs(All_Start(i,3)-Temp(3))+abs(All_Start(i,4)-Temp(4)))/4<10
            Temp=[];
            break;
        end
    end
    All_Start=[All_Start;Temp];
end
end

