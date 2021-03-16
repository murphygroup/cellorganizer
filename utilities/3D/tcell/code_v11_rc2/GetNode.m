function [JiaoDian_Point,Idexk2]=GetNode(SmallCellData,radii,centers)
syms x y
JiaoDian_Point=[];
%  if the small cells and big cells have the intersection ,it will be saved and outputed 
Idexk2=[];
centers=round(centers);
for k2=1:length(radii)
    % if the center of the ellipse and the center distance of the large circle are larger than the sum of the ellipse's long axis and the radius of the large circle.
    % there is no intersection between the two circles
    tmp=sqrt(   (centers(k2,1)-SmallCellData.Center(1))^2+(centers(k2,2)-SmallCellData.Center(2))^2  );
    if (tmp<(max(abs(SmallCellData.Zhou))/2+radii(k2)))&&(tmp>radii(k2))
        p=SmallCellData.equ{2};
        [x_solve,y_solve]=vpasolve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0, (x-centers(k2,1))^2+(y-centers(k2,2))^2== radii(k2)^2, x, y, [-Inf Inf]);
        solutions2=double([x_solve,y_solve]);
        solutions2 = solutions2(imag(solutions2(:, 1)) == 0, :);
        if ~isempty(solutions2)
                JiaoDian_Point=[JiaoDian_Point;solutions2];
                Idexk2=[Idexk2;k2];
        end
    end
end