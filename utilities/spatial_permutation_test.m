function p_val = spatial_permutation_test(obs, x, y, z)
d = length(x);

d1 = fix(d/2);
d2 = d - d1;
score_list = [];

for i=1:50
    idx1 = randperm(d);
    x1 = x(idx1);
    y1 = y(idx1);
    z1 = z(idx1);
  
    [X,Y,Z] = sphere;
    [az,phi,r]=cart2sph(X(:),Y(:),Z(:));
    new_az=az;
    new_phi=phi;
    new_r=r;

    az_ = [];
    phi_ = [];
    r_ = [];
    for i=1:10
        az_ = [az_ az'];
        phi_ = [phi_ phi'];
        r_ = [r_ 0.1*i*r'];
    end
    [new_x,new_y,new_z]=sph2cart(az_',phi_',r_');

    y_b1 = mvksdensity([x1(1:d1)' y1(1:d1)' z1(1:d1)'],[new_x new_y new_z],'Bandwidth',0.28,'Kernel','normal');
    y_b2 = mvksdensity([x1(d1+1:end)' y1(d1+1:end)' z1(d1+1:end)'],[new_x new_y new_z],'Bandwidth',0.28,'Kernel','normal');


    y_b1 = y_b1 + 1e-5;
    y_b2 = y_b2 + 1e-5;
    y_b1 = y_b1/sum(y_b1);
    y_b2 = y_b2/sum(y_b2);
    KL_div = sum(y_b1.*log(y_b1./y_b2));
    score_list(end+1) = KL_div;
end

p_val = sum(obs<score_list)/50;
end