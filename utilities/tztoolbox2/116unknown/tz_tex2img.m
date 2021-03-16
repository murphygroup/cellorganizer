function img=tz_tex2img(tex,templ,imgsize)

%function img=tz_tex2img(tex,templ,imgsize)

if isempty(templ)
    templ=zeros(imgsize);
end

ptex=tex/sum(tex(:));

diagtex=diag(diag(ptex));
triutex=triu(ptex);
symtex=triutex+triutex'-diagtex;

imgsize=size(templ);
img=zeros(imgsize);

mtex=sum(tex+tex',1)/sum(tex(:))/2;
tmtex=mtex;
tmtex(mtex==0)=1;
ctex=symtex./repmat(tmtex,size(tex,1),1);

ptex=symtex;
%four corners
neighbors=[templ(1,2),templ(2,1)];
neighbors(neighbors==0)=0;
p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
img(1,1)=p;
templ(1,1)=p;

neighbors=[templ(1,end-1),templ(2,end)];
neighbors(neighbors==0)=0;
p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
img(1,end)=p;
templ(1,end)=p;

neighbors=[templ(end,2),templ(end-1,1)];
neighbors(neighbors==0)=0;
p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
img(end,1)=p;
templ(end,1)=p;

neighbors=[templ(end,end-1),templ(end-1,end)];
neighbors(neighbors==0)=0;
p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
img(end,end)=p;
templ(end,end)=p;

%four borders
for i=2:imgsize(2)-1
    neighbors=[templ(1,i-1),templ(1,i+1),templ(2,i)];
    neighbors(neighbors==0)=0; 
    p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
    img(1,i)=p;
    templ(1,i)=p;
end

for i=2:imgsize(1)-1
    neighbors=[templ(i-1,1),templ(i+1,1),templ(i,2)];
    neighbors(neighbors==0)=0; 
    p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
    img(i,1)=p;
    templ(i,1)=p;
end

for i=2:imgsize(2)-1
    neighbors=[templ(end,i-1),templ(end,i+1),templ(end-1,i)];
    neighbors(neighbors==0)=0; 
    p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
    img(end,i)=p;
    templ(end,i)=p;
end

for i=2:imgsize(1)-1
    neighbors=[templ(i-1,end),templ(i+1,end),templ(i,end-1)];
    neighbors(neighbors==0)=0; 
    p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
    img(i,end)=p;
    templ(i,end)=p;
end

for i=2:imgsize(1)-1
    for j=2:imgsize(2)-1
        neighbors=[templ(i-1,j-1),templ(i+1,j+1),templ(i-1,j+1),templ(i+1,j-1)];
        neighbors(neighbors==0)=0; 
        p=tz_tex2pixel(neighbors,ptex,mtex,ctex);
        
        img(i,j)=p;
        templ(i,j)=p;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=tz_tex2pixel(neighbors,ptex,mtex,ctex)

%function p=tz_tex2pixel(neighbors,tex)

sr=neighbors;
sr(neighbors==0)=[];
if isempty(sr)
    p=find(tz_mnornd(1,mtex(1:end-1),1));
else
%     tmtex=mtex;
%     tmtex(mtex==0)=1;
    bterm=ptex(:,sr(1))';
    for i=1:length(sr)-1
        bterm=bterm.*ptex(:,sr(i))';
    end
    
    comtex=bterm;
    
    if(sum(comtex)==0)
        p=0;
    else
        comtex=comtex/sum(comtex);
        
        p=find(tz_mnornd(1,comtex(1:end-1),1));
    end

end
