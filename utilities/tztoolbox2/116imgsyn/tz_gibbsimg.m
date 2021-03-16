function img = tz_gibbsimg(imgsize,niter)
%TZ_GIBBSIMG Under construction.

%function img = tz_gibbsimg(imgsize)
%OVERVIEW
%   Synthesize an image by Gibbs sampling
%PARAMETERS
%   imgsize - the size of the output image
%   niter - number of iteration
%RETURN
%   img - synthesized image
%DESCRIPTION
%   
%HISTORY
%   29-Mar-2005 Initial write TINGZ
%SEE ALSO
%   

warning('Under developing ...');
p=0.5;

img=binornd(255,p,imgsize(1)+2,imgsize(2)+2);

for k=1:niter
    k
    for pt=0:1
        for i=2:imgsize(1)+1
            s=mod(i+pt,2)+2;
            for j=s:2:imgsize(2)+1
                nbs=[img(i-1,j),img(i,j-1),img(i+1,j),img(i,j+1)];
                ws=[0.1 0.2 0.1 0.2]/100;
                img(i,j)=binornd(255,mean(nbs)/255,1,1);
            end
        end
    end
%     for i=2:imgsize(1)+1
%         s=mod(i+1,2)+2;
%         for j=s:2:imgsize(2)
%             nbs=[img(i-1,j),img(i,j-1),img(i+1,j),img(j+1)];
%             img(i,j)=binornd(255,exp(-sum(nbs.*[0.1 0.2 0.1 0.2])),1,1);
%         end
%     end
end
