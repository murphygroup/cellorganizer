% 2011-11-18 tebuck: adding 'nucStats' set.
% 2011-12-01 tebuck: adding 'protTotalIntensity' set.

% $$$ fprintf('##########\nComputing set ''%s''\n\n', fsetnames{zed})
% $$$ fprintf('########## Computing set ''%s''\n', fsetnames{zed})

switch fsetnames{zed}
  case 'tas'
% $$$     [names, feats, slfnames] = ml_tas( protstruct.channel,0);

% $$$     if optimize && any(strcmpi(naming_convention.blank_channels, 'protein'))
% $$$       [names, feats, slfnames] = ml_tas(0, 0);
% $$$       feats = nan(size(feats)); 
% $$$     else
% $$$       [names, feats, slfnames] = ml_tas(, 0);
% $$$     end
   
    prot = protstruct.channel;
    commonScriptTAS
  
  case 'tasx2'
% $$$     downsampleProteinChannel
% $$$     
% $$$     [names, feats, slfnames] = ml_tas( protstruct.downsampled2x,0);
    
% $$$     if any(strcmpi(naming_convention.blank_channels, 'protein'))
% $$$       [names, feats, slfnames] = ml_tas(0, 0);
% $$$       feats = nan(size(feats)); 
% $$$     else
% $$$       downsampleProteinChannel
% $$$       [names, feats, slfnames] = ml_tas(protstruct.downsampled2x, 0);
% $$$     end

    downsampleProteinChannel
   
    prot = protstruct.downsampled2x;
    commonScriptTAS
   
    renameDownsampled

  case 'mutualInfo'
    loadNucleusChannel; loadMicrotubuleChannel; loadERChannel; 
    
% $$$     sprot = size(protstruct.channel);
% $$$     prot_re = reshape(protstruct.channel,[sprot(1)*sprot(2) 1]);
% $$$     nuc_re = reshape(nucstruct.channel,[sprot(1)*sprot(2) 1]);
% $$$     tub_re = reshape(tubstruct.channel,[sprot(1)*sprot(2) 1]);
% $$$     er_re = reshape(erstruct.channel,[sprot(1)*sprot(2) 1]);
    prot_re = reshape(protstruct.channel, [], 1);
    nuc_re = reshape(nucstruct.channel, [], 1);
    tub_re = reshape(tubstruct.channel, [], 1);
    er_re = reshape(erstruct.channel, [], 1);

    commonScriptMutualInformation
    
  case 'mutualInfox2'
    loadNucleusChannel; loadMicrotubuleChannel; loadERChannel; 

    downsampleAllChannels

% $$$     sprot = size(protstruct.downsampled2x);
% $$$     prot_re = reshape(protstruct.downsampled2x,[sprot(1)*sprot(2) 1]);
% $$$     nuc_re = reshape(nucstruct.downsampled2x,[sprot(1)*sprot(2) 1]);
% $$$     tub_re = reshape(tubstruct.downsampled2x,[sprot(1)*sprot(2) 1]);
% $$$     er_re = reshape(erstruct.downsampled2x,[sprot(1)*sprot(2) 1]);
    prot_re = reshape(protstruct.downsampled2x, [], 1);
    nuc_re = reshape(nucstruct.downsampled2x, [], 1);
    tub_re = reshape(tubstruct.downsampled2x, [], 1);
    er_re = reshape(erstruct.downsampled2x, [], 1);

    commonScriptMutualInformation

    renameDownsampled

 case 'texture'
% $$$     fprintf('#####\nIn case ''texture''\n\n')
% $$$     protstruct
    
    prot = protstruct.channel;
    commonScriptTexture

  case 'texturex2'
% $$$     fprintf('#####\nIn case ''texturex2''\n\n')
% $$$     protstruct
    
    downsampleProteinChannel
% $$$     protstruct

    prot = protstruct.downsampled2x;
    commonScriptTexture

    renameDownsampled

  case 'texturex4'
% $$$     fprintf('#####\nIn case ''texturex4''\n\n')
% $$$     protstruct
    
    DSF = 4;
    downsampleProteinChannel
% $$$     protstruct

    prot = protstruct.downsampled4x;
    commonScriptTexture

    renameDownsampled
    DSF = 2;

  case 'overlap'
    loadNucleusChannel; loadMicrotubuleChannel; loadERChannel;

    thresholdAllChannels

    tmp_prot = protstruct.channel;
    tmp_prot_fg = protstruct.channel_fg;
    tmp_nuc_fg = nucstruct.channel_fg;
    tmp_tub_fg = tubstruct.channel_fg;
    tmp_er_fg = erstruct.channel_fg;
    tmp_prot_mfg = protstruct.channel_mfg;
    tmp_prot_large_fg = protstruct.channel_large_fg;
    tmp_prot_large_mfg = protstruct.channel_large_mfg;

    commonScriptOverlap

  case 'overlapx2'
    loadNucleusChannel; loadMicrotubuleChannel; loadERChannel;

    downsampleAllChannels

    thresholdDownsampled2xAllChannels

    tmp_prot = protstruct.downsampled2x;
    tmp_prot_fg = protstruct.downsampled2x_fg;
    tmp_nuc_fg = nucstruct.downsampled2x_fg;
    tmp_tub_fg = tubstruct.downsampled2x_fg;
    tmp_er_fg = erstruct.downsampled2x_fg;
    tmp_prot_mfg = protstruct.downsampled2x_mfg;
    tmp_prot_large_fg = protstruct.downsampled2x_large_fg;
    tmp_prot_large_mfg = protstruct.downsampled2x_large_mfg;

    commonScriptOverlap

    renameDownsampled

 case 'nonObjFluor'
    if optimize
      thresholdProteinChannel
    else
      loadNucleusChannel; loadMicrotubuleChannel; loadERChannel;

      thresholdAllChannels
    end
  
    prot_int = sum(protstruct.channel_fg(:));
    bg_prot = protstruct.channel;
    bg_prot(protstruct.channel>protstruct.channel_thr) = 0;
    prot_int_large = sum(protstruct.channel_large_fg(:));

    prot_int_mthr = sum(protstruct.channel_mfg(:));
    bg_prot_mthr = protstruct.channel;
    bg_prot_mthr(protstruct.channel>protstruct.channel_mthr) = 0;
    prot_int_large_mthr = sum(protstruct.channel_large_mfg(:));

    commonScriptNOF

  case 'nonObjFluorx2'
    if optimize
      downsampleProteinChannel
      thresholdDownsampled2xProteinChannel
    else
      loadNucleusChannel; loadMicrotubuleChannel; loadERChannel;

      downsampleAllChannels

      thresholdDownsampled2xAllChannels
    end

    prot_int = sum(protstruct.downsampled2x_fg(:));
    bg_prot = protstruct.downsampled2x;
    bg_prot(protstruct.downsampled2x>protstruct.downsampled2x_thr) = 0;
    prot_int_large = sum(protstruct.downsampled2x_large_fg(:));

    prot_int_mthr = sum(protstruct.downsampled2x_mfg(:));
    bg_prot_mthr = protstruct.downsampled2x;
    bg_prot_mthr(protstruct.downsampled2x>protstruct.downsampled2x_mthr) = 0;
    prot_int_large_mthr = sum(protstruct.downsampled2x_large_mfg(:));

    commonScriptNOF

    renameDownsampled

  case 'obj'
    loadNucleusChannel

    thresholdProteinChannel;  thresholdNucleusChannel;

    protobjs = protstruct.channel_objectsizes;
    largeprotobjs = protstruct.channel_large_objectsizes;

    protobjs_mthr = protstruct.channel_mobjectsizes;
    largeprotobjs_mthr = protstruct.channel_large_mobjectsizes;

    nucobjs = nucstruct.channel_objectsizes;

    normEuN = bweuler(protstruct.channel_fg>0,8)/length(nucobjs);
    normEuN_mthr = bweuler(protstruct.channel_mfg>0,8)/length(nucobjs);

    commonScriptObject;

  case 'objx2'
    loadNucleusChannel

    downsampleProteinChannel;  downsampleNucleusChannel;

    thresholdDownsampled2xProteinChannel;  thresholdDownsampled2xNucleusChannel;

    protobjs = protstruct.downsampled2x_objectsizes;
    largeprotobjs = protstruct.downsampled2x_large_objectsizes;

    protobjs_mthr = protstruct.downsampled2x_mobjectsizes;
    largeprotobjs_mthr = protstruct.downsampled2x_large_mobjectsizes;

    nucobjs = nucstruct.downsampled2x_objectsizes;

    normEuN = bweuler(protstruct.downsampled2x_fg>0,8)/length(nucobjs);
    normEuN_mthr = bweuler(protstruct.downsampled2x_mfg>0,8)/length(nucobjs);

    commonScriptObject;

    renameDownsampled

  case 'objRegion'
    loadNucleusChannel

    thresholdProteinChannel;  thresholdNucleusChannel;

% $$$     [names, feats, slfnames] = ml_imgfeatures( protstruct.channel_fg, nucstruct.channel_fg);
% $$$ 
% $$$     if length(feats)~=14
% $$$         feats = repmat(nan, [1 14]);
% $$$     end
    
    protein_image = protstruct.channel_fg;
    nuclear_image = nucstruct.channel_fg;
    
    commonScriptObjectRegion

  case 'objRegionx2'
    loadNucleusChannel

    downsampleProteinChannel;  downsampleNucleusChannel;

    thresholdDownsampled2xProteinChannel;  thresholdDownsampled2xNucleusChannel;

% $$$     [names, feats, slfnames] = ml_imgfeatures( protstruct.downsampled2x_fg, nucstruct.downsampled2x_fg);
% $$$ 
% $$$     if length(feats)~=14
% $$$         feats = repmat(nan, [1 14]);
% $$$     end

    protein_image = protstruct.downsampled2x_fg;
    nuclear_image = nucstruct.downsampled2x_fg;
    
    commonScriptObjectRegion

    renameDownsampled


  case 'nucStats'
    loadNucleusChannel

    thresholdNucleusChannel;

    nucobjs = nucstruct.channel_objectsizes;

    commonScriptNucStats

    
  case 'protTotalIntensity'
    %loadNucleusChannel

    commonScriptProtTotalIntensity

    
  case 'protMeanIntensity'
    loadMicrotubuleChannel;
    
    thresholdAllChannels;
    
    commonScriptProtMeanIntensity

    
  otherwise
    error('improper feature set specified');
end
