function consolidated_model = consolidateshapespacemodels(model1,model2)
% input models should be the cellShapeModels substruct

X1=model1.X;
coeff1=model1.coeff;
score1=model1.train_score; 

X2=model2.X;
coeff2=model2.coeff;
score2=model2.train_score; 

consolidateX=[X1;X2];

%default value in cellorganizer
latent_dim=100;

[scales,coeff,score,latent,tsquared,explained,mu,train_score_consolidate,train_explained,train_coeff_consolidate] = CalculatePCA(consolidateX, latent_dim);
spharm1=model1.all_spharm_descriptors;
spharm2=model2.all_spharm_descriptors;

%Create new model with required features 
consolidated_model.coeff=train_coeff_consolidate;
consolidated_model.X=consolidateX;
consolidated_model.train_score=train_score_consolidate;
consolidated_model.numimgs1=model1.numimgs;
consolidated_model.numimgs2=model2.numimgs;
consolidated_model.numimgs=model1.numimgs+model2.numimgs;
consolidated_model.all_spharm_descriptors=cat(3,spharm1,spharm2);
consolidated_model.hausdorff_distances=horzcat(model1.hausdorff_distances,model2.hausdorff_distances);
if ismember(model1.components,model2.components)
    consolidated_model.components=model2.components;
else
    consolidated_model.components={};
    consolidated_model.components={model1.components,model2.components};
end
