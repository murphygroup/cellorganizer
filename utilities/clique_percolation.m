function [ KLD ] = clique_percolation( models)
%Author : Serena Abraham

spatial1=models{1};
spatial2=models{2};

graphSize = 50;

for a=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/graphSize)
    for b=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/graphSize)
        normdists=horzcat(spatial1.normdists((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.normdists((a-1)*graphSize+1:(a-1)*graphSize+graphSize));
        anglestheta=horzcat(spatial1.anglestheta((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.anglestheta((a-1)*graphSize+1:(a-1)*graphSize+graphSize));
        anglesphi=horzcat(spatial1.anglesphi((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.anglesphi((a-1)*graphSize+1:(a-1)*graphSize+graphSize));

        [x,y,z]=sph2cart(anglestheta,anglesphi,normdists);
        calc_avg=[];

        distances=zeros(length(x),length(x));
        for i=1:length(x)
            for j=1:length(x)
                if distances(i,j)==0
                    if i==j
                        dist=0;
                        distances(i,j)=dist;
                    else
                        dist=norm([x(i),y(i),z(i)]-[x(j),y(j),z(j)]);
                        distances(i,j)= dist;
                        distances(j,i)=dist;   
                    end
                    calc_avg(end+1)=dist;
                end
            end
        end
    end
end
thresh=mean(calc_avg);
prob1=[];
prob2=[];

% idx1 = randperm(length(spatial1.normdists));
% spatial1.normdists = spatial1.normdists(idx1);
% spatial1.anglestheta = spatial1.anglestheta(idx1);
% spatial1.anglesphi = spatial1.anglesphi(idx1);
% idx2 = randperm(length(spatial2.normdists));
% spatial2.normdists = spatial2.normdists(idx2);
% spatial2.anglestheta = spatial2.anglestheta(idx2);
% spatial2.anglesphi = spatial2.anglesphi(idx2);

for a=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/graphSize)
    %for b=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
        normdists=horzcat(spatial1.normdists((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.normdists((a-1)*graphSize+1:(a-1)*graphSize+graphSize));
        anglestheta=horzcat(spatial1.anglestheta((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.anglestheta((a-1)*graphSize+1:(a-1)*graphSize+graphSize));
        anglesphi=horzcat(spatial1.anglesphi((b-1)*graphSize+1:(b-1)*graphSize+graphSize),spatial2.anglesphi((a-1)*graphSize+1:(a-1)*graphSize+graphSize));

        [x,y,z]=sph2cart(anglestheta,anglesphi,normdists);

        distances=zeros(length(x),length(x));
        for i=1:length(x)
            for j=1:length(x)
                if distances(i,j)==0
                    if i==j
                        dist=0;
                        distances(i,j)=dist;
                    else
                        dist=norm([x(i),y(i),z(i)]-[x(j),y(j),z(j)]);
                        distances(i,j)= dist;
                        distances(j,i)=dist;   
                    end
                end
            end
        end
        
        
        graph=zeros(length(x),length(x));
        for i=1:length(x)
            for j=1:length(x)
                if distances(i,j)>=thresh
                    graph(i,j)=1;
                end
            end
        end
       c2=ELSclique(graph);
       %Traverse cliques
       count1=0;
       count2=0;
       c2 = full(c2);
       [~, d2] = size(c2);
       for k=1:d2
           clique=find(c2(:, k));
           if length(clique) > 1
                for l=1:length(clique)
                    if clique(l)<= v
                        count1=count1+1;
                    else
                        count2=count2+1;
                    end
                end
            prob1(end+1)=count1/(count1+count2);
            prob2(end+1)=count2/(count1+count2);
            count1=0;
            count2=0;
            end
        end
    %end
end

%KL Divergence
% KLD1=KLDiv(prob1',prob2');
prob1 = prob1 + 1e-7;
prob2 = ones(size(prob1))/2;
% prob1 = prob1/sum(prob1);
% prob2 = prob2/sum(prob2);
KLD = mean(prob1.*log(prob1./prob2) + (1-prob1).*log((1-prob1)./(1-prob2)))
end




