function W = makecube(imgsize, s);
[aa, bb, cc] = meshgrid(1:imgsize, 1:imgsize, 1:imgsize);
space = (imgsize - s)/2;
Q = aa <= imgsize-space;
R = aa > space;

A = Q & R;

S = bb <= imgsize-space;
T = bb > space;

B = S & T;

U = cc <= imgsize-space;
V = cc > space;

C = U & V;

W = A & B & C;
end%makecube
