function [objgauss,volume] = tp_gengaussdistr(model,n)
% Generate a set of gaussian distributions from the model

% Graphical Model: Y <-- X --> Z

objgauss = cell(n,1);
i = 0;
while i < n
    % Generate covariance matrices
    x = ml_rnd(model.x);
    model.y_x.mu = model.y_x.a1*(1-exp(-model.y_x.b1*x));
    model.z_x.mu = model.z_x.a1*(1-exp(-model.z_x.b1*x));
    model.y_x.sigma = model.y_x.a2*(1-exp(-model.y_x.b2*x));
    model.z_x.sigma = model.z_x.a2*(1-exp(-model.z_x.b2*x));
    y = ml_rnd(model.y_x);
    z = ml_rnd(model.z_x);

    if y>x || z>x || z>y || x<=0 || y<=0 || z<=0
        continue
    end
    i = i+1;
%     newSigma = diag([x y z].^2)/25;
    newSigma = diag([x y z].^2);    
    volume(i) = x*y*z;

    % Random rotation
    t = rand * pi;
    Rx = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
    t = rand * pi;
    Ry = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
    t = rand * pi;
    Rz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
    R = Rz*Ry*Rx;
    newCov = R*newSigma*R';
    objgauss{i} = newCov;
end
