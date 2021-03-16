function inst_shape_model

%% Calculate distance metric

%% MDS shape space dimension reduction
[Y,R] = mdsClassic(D,1:5);
X = Y.coords{end}';

%% Bandwidth selection
h = opt_bandwd(X);

%% Triangulate space
tes = delaunayn(X);

%% Generate sample point
xx = kern_rnd(X,h);

%% Locate sample point in triangle mesh
simplex_idx = tsearchn(X,tes,xx);
xidx = tes(simplex_idx,:);
xtri = X(xidx,:);

%% Shape interpolation path
[lambda,ordering] = shape_interp_track(xtri',xx');
interporder = fliplr(xidx(ordering));
interplambda = fliplr(lambda);