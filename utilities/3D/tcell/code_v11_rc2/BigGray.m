function [ Avg_Gray ] = BigGray( Superposition,k,b,TempLeft,TempRight )
[m, n] = size(Superposition);
if any(TempLeft >= [n, m] + 0.5) || any(TempLeft <= -0.5) || any(TempRight >= [n, m] + 0.5) || any(TempRight <= -0.5)
    Avg_Gray = 0;
    return;
end
    
DisX=abs(TempRight(1)-TempLeft(1));
DisY=abs(TempRight(2)-TempLeft(2));
if(DisX>DisY)
    if(TempLeft(1)>TempRight(1))
        Temp=TempRight;
        TempRight=TempLeft;
        TempLeft=Temp;
    end
    x=TempLeft(1):0.01:TempRight(1);
    y=k*x+b;
    Gray=0;
    X=round(x);
    Y=round(y);
        
    nnn=find(X>n | X<1);
    X(nnn)=[];
    Y(nnn)=[];    
    mmm=find(Y>m | Y<1);
    X(mmm)=[];
    Y(mmm)=[];  
    
    for i=1:length(X)
        Gray=Gray+double(Superposition(Y(i),X(i)));
    end
    Avg_Gray=Gray/length(X);
else
    if(TempLeft(2)>TempRight(2))
        Temp=TempRight;
        TempRight=TempLeft;
        TempLeft=Temp;
    end
    y=TempLeft(2):0.01:TempRight(2);
    x=(y-b)/k;
    Gray=0;
    X=round(x);
    Y=round(y);
    
    nnn=find(X>n | X<1);
    X(nnn)=[];
    Y(nnn)=[];    
    mmm=find(Y>m | Y<1);
    X(mmm)=[];
    Y(mmm)=[];  
    
    for i=1:length(X)
        Gray=Gray+double(Superposition(Y(i),X(i)));
    end
    Avg_Gray=Gray/length(X);
end

