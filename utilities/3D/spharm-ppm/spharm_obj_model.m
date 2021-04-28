function model = spharm_obj_model(options)
% 1/26/2021 R.F. Murphy - save normalized nuclear distances into model
% (only for objects that were successfully parameterized)
% 3/2/2021 R.F. Murphy - remove unused ppm code
% 4/24/2021 R.F. Murphy - merge in Serena's changes

    options_spharm = options.options_spharm;
    spharm_obj_dir = [pwd filesep 'spharm_input'];
    param_save_dir = [pwd filesep 'param'];
    spharm_obj_files = ml_ls(spharm_obj_dir);
    for i=1:length(spharm_obj_files)
        spharm_obj_files{i}=[spharm_obj_dir filesep spharm_obj_files{i}];
    end
    disp('build spharm model for spharm object');
    [status, msg, msgID] = mkdir('spharm_result');
    cd spharm_result
    tic; %answer = img2slml2('3D', spharm_obj_files, spharm_obj_files, [], options_spharm); 
    answer = spharm_rpdm_v2(spharm_obj_files,options_spharm);
    toc,
    
    t = load([options_spharm.model.filename(1:end-3) 'mat']);
    t.model.cellShapeModel.numimgs
    t.model.spatial = get_obj_coordinates(param_save_dir, ...
        t.model.cellShapeModel.parameterization_successful);
    model.spharm_obj_model = t.model;
    cd ..
end

function spatial = get_obj_coordinates(paramdir,isdone)

spatial = [];
normdists = [];
anglestheta = [];
anglesphi = [];
mappos_x = [];
distcodes=[];
paramfiles = ml_ls([paramdir filesep '*.mat']);
for i=1:length(paramfiles)
    pp = load(paramfiles{i});
%    for j length(prot.seg.centers)
%        all_indices = [all_indices; [i,j]];
    normdists = [normdists; pp.prot.seg.normdists];
    anglestheta = [anglestheta; pp.prot.seg.angles.theta];
    anglesphi = [anglesphi; pp.prot.seg.angles.phi];
    mappos_x = [mappos_x; pp.prot.seg.mappos_x];
    distcodes=[distcodes;pp.prot.seg.distcodes];
end
spatial.normdists = normdists(isdone)';
spatial.anglestheta = anglestheta(isdone)';
spatial.anglesphi = anglesphi(isdone)';
spatial.mappos_x = mappos_x(isdone,:);
spatial.distcodes=distcodes(isdone,:);
beta=ml_logreg(spatial.mappos_x(:,2:6),spatial.distcodes(:,3));
spatial.beta=beta;
end