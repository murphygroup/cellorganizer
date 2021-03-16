function cellShapeModel = param2cell_ratio_model_2D(cellparam, cellshape_options)

if length( cellparam ) >= 10
    number_of_components = 10;
else
    number_of_components = length( cellparam );
end


cellshape_options = ml_initparam( cellshape_options, ...
    struct('modelname','rds', ...
    'screencell','no',...
    'ml_rdistpca',struct('ml_estpdf', ...
    struct('name','mvn','transform',struct('funname','_pca', ...
    'param',struct('ncomp', number_of_components))), ...
    'startangle','nuc')));

%hidden parameter

cellshape_options.anglestep = 1;

cellIndices = 1:length(cellparam);



combfeats = ml_combccfeats(cellparam(cellIndices));

if strcmp(cellshape_options.screencell,'yes')==1
    %train cell shape
    %Area ratios
    % combratio = combfeats(goodCellIndices,2);
    combratio = combfeats(:,2);
    
    %Gaussian mixture of area ratios
    [bestk,bestpp,bestmu,bestcov] = ...
        ml_mixtures4(combratio',1,3,0,1e-5,0);
    
    if bestk>1
        %         pp1=bestpp(1);
        %         pp2=bestpp(2);
        sigma1=sqrt(bestcov(:,:,1));
        sigma2=sqrt(bestcov(:,:,2));
        mu1=bestmu(1);
        mu2=bestmu(2);
        
        %Find threshold for small and big ratios
        x = ml_gscross(mu1,sigma1,mu2,sigma2);
        
        sizeThreshold = x(1);
        
        smallCellIndices = cellIndices(combratio<sizeThreshold);
        %         largeCellIndices = cellIndices(combratio>=sizeThreshold);
    else
        smallCellIndices = cellIndices;
    end
else
    smallCellIndices = cellIndices;
end

%PCA ratio model
%model.cellShapeModel = tz_rdistpca(cellcodes(smallCellIndices),param);

cellparam = cellparam(smallCellIndices);

%Cell shape model
%icaoberg april 11, 2012
cellShapeModel = ml_trainshapemodel2D(cellparam, ...
    cellshape_options );

