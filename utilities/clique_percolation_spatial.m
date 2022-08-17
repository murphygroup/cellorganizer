function [ KLD ] = clique_percolation_spatial(x1, x2)
%Author : Serena Abraham

graphSize = 50;

for a=1:floor(min(length(x1),length(x2))/graphSize)
    for b=1:floor(min(length(x1),length(x2))/graphSize)
        temp=vertcat(x1((b-1)*graphSize+1:(b-1)*graphSize+graphSize, :),...
            x2((a-1)*graphSize+1:(a-1)*graphSize+graphSize, :));
        

        x = temp(:, 1);
        y = temp(:, 2);
        z = temp(:, 3);
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

% idx1 = randperm(length(x1));
% x1 = x1(idx1, :);
% 
% idx2 = randperm(length(x2));
% x2 = x2(idx2, :);

for a=1:floor(min(length(x1),length(x2))/graphSize)
    %for b=1:floor(min(length(spatial1.normdists),length(spatial2.normdists))/20)
        temp=vertcat(x1((a-1)*graphSize+1:(a-1)*graphSize+graphSize, :),...
            x2((a-1)*graphSize+1:(a-1)*graphSize+graphSize, :));
        
        x = temp(:, 1);
        y = temp(:, 2);
        z = temp(:, 3);

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
        
        
        graph_=zeros(length(x),length(x));
        for i=1:length(x)
            for j=1:length(x)
                if distances(i,j)>=thresh
                    graph_(i,j)=1;
                end
            end
        end
       if sum(graph_, 'all') == 0
           continue
       end
       c2=ELSclique(graph_);

       %Traverse cliques
       count1=0;
       count2=0;
       c2 = full(c2);
       [~, d2] = size(c2);
       for k=1:d2
           clique=find(c2(:, k));
           if length(clique) > 2
                for l=1:length(clique)
                    if clique(l)<=graphSize
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
KLD = mean(prob1.*log(prob1./prob2) + (1-prob1).*log((1-prob1+ 1e-7)./prob2))
end




