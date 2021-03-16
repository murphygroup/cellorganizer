function SmallCellData=Process_SmallCell(img,loc,width,height)
% width��height�ǿ�ȡС���ο�Ŀ�Ⱥ͸߶�
position=[loc(1)-width/2,loc(2)-height/2,width,height];
% �߽紦��  500��444�ֱ���ͼ��Ŀ�͸�
if position(1)<1
    position(1)=1;
end
if position(1)+width>size(img,2)
    position(1)=size(img,2)-width;
end
if position(2)<1
    position(2)=1;
end
if position(2)+height>size(img,1)
    position(2)=size(img,1)-height;
end
tmp=imcrop(img,position);
BW=im2bw(tmp,graythresh(tmp));
% �õ������ͨ�� ������Ϊ0
L = bwlabel(BW);%�����ͨ����
stats = regionprops(L);
Ar = cat(1, stats.Area);
ind = find(Ar ==max(Ar));%�ҵ������ͨ����ı��
if(length(ind)==1)
    BW(L~=ind)=0;%������������Ϊ0
    MAX_AR=max(Ar);
else
    MAX_AR=0;
end

% ����ȡ�������е������ͨ�����ռ���ο����ı�ֵ����0.5 ˵�����Ŀ�����Ч
if MAX_AR/width*height>0.3
    BW=medfilt2(BW,[2 2]);
    % Edge=edge(BW,'canny');
    % for bw image, use bwboundaries to get boundary points
    BW = imfill(BW, 'hole');
    B = bwboundaries(BW);
    if numel(B) == 1
        m = B{1}(:, 1);
        n = B{1}(:, 2);
    else
        SmallCellData = struct();
        return;
    end
    % % ��ȡ�߽� ��Բ���
    % [m,n]=find(Edge==1);
    % % ����껹ԭ��ԭͼ��
    m=m+loc(2)-height/2;n=n+loc(1)-width/2;
    SmallCellData=Ellipse_Fitting2([n,m]);
else
    SmallCellData=struct();
end

    