function [ output_args ] = clique_percolation( models,fileID, options)
%Author : Serena Abraham

spatial1=models{1};
spatial2=models{2};

for a=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
    for b=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
        normdists=horzcat(spatial1.normdists((b-1)*20+1:(b-1)*20+20),spatial2.normdists((a-1)*20+1:(a-1)*20+20));
        anglestheta=horzcat(spatial1.anglestheta((b-1)*20+1:(b-1)*20+20),spatial2.anglestheta((a-1)*20+1:(a-1)*20+20));
        anglesphi=horzcat(spatial1.anglesphi((b-1)*20+1:(b-1)*20+20),spatial2.anglesphi((a-1)*20+1:(a-1)*20+20));

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

for a=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
    %for b=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
        normdists=horzcat(spatial1.normdists((b-1)*20+1:(b-1)*20+20),spatial2.normdists((a-1)*20+1:(a-1)*20+20));
        anglestheta=horzcat(spatial1.anglestheta((b-1)*20+1:(b-1)*20+20),spatial2.anglestheta((a-1)*20+1:(a-1)*20+20));
        anglesphi=horzcat(spatial1.anglesphi((b-1)*20+1:(b-1)*20+20),spatial2.anglesphi((a-1)*20+1:(a-1)*20+20));

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
       [c1,c2,c3]=k_clique(3,graph);
       %Traverse cliques
       count1=0;
       count2=0;
       for k=1:length(c2)
           clique=c2{k};
            for l=1:length(clique)
                if clique(l)<=20
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
    %end
end

%KL Divergence
P1=unique(prob1)';
P1=P1/sum(P1);
P2=unique(prob1)';
P2=P2/sum(P2);
X=(1:length(P2))';
output_args=KLDiv(X,P1,P2);




