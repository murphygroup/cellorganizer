function [mse_std, sse, pos2_align ] = procrustes_mse( pos1, pos2 )
%function for determining the MSE between points. The points pos2 are alligned to pos1 via procrustes superimposition and the MSE is taken 

d1 = size(pos1,2);
d2 = size(pos2,2);

if d1-d2 < 0
    pos1 = padarray(pos1, [0,abs(d1-d2)], 'post');
else
    pos2 = padarray(pos2, [0, d1-d2], 'post');
end
    
[mse_std, pos2_align] = procrustes(pos1, pos2);    

sse = sum((pos1(:) - pos2_align(:)).^2);


end

