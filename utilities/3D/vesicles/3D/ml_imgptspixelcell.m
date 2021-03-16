function [point,distance] = ml_imgptspixelcell( cellslice,pts,nuccenter )

ptstrim=pts;
imgsize=size(cellslice);
ptstrim(ptstrim(:,1)<=0 | ptstrim(:,1)>imgsize(1),:)=[];
ptstrim(ptstrim(:,2)<=0 | ptstrim(:,2)>imgsize(2),:)=[];

if isempty(ptstrim)
    error('empty points');
end

%imshow(cellslice); hold on;
%plot(ptstrim(:,1),ptstrim(:,2),'r-');

idx=sub2ind(imgsize,ptstrim(:,1),ptstrim(:,2));

ps=cellslice(idx);
intc2=find(ps>0);

if ~isempty(intc2)
%    plot(ptstrim(intc2,1),ptstrim(intc2,2),'g+'); hold off; pause(0.1);
            
    distances = sqrt((ptstrim(intc2,1)-nuccenter(1)).^2 + ... 
            (ptstrim(intc2,2)-nuccenter(2)).^2);
    [distance,minidx]=min(distances);
    point=ptstrim(intc2(minidx),:);
else
    warning('ml_imgptspixelcell: No intercept.')
    [point, distance] = findClosest( cellslice, pts );
end

end

