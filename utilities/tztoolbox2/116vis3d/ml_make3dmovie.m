function mov = ml_make3dmovie(h)

if ~exist('h','var')
    h = gcf;
end

axis tight
set(gca,'CameraViewAngleMode','manual')
viewAngles = 0:10:350;

for k=1:length(viewAngles)
    view(10,[viewAngles(k)]);
    mov(k) = getframe(h,[1 1 436 108*4]);
end

