function [] = POSTGenHybridComp( LINKSDat, fibre1Dat, fibre2Dat, matrixDat,...
    RealizationFolder, MeshFolder, s_info )


%% Read LINKS input options
[~,s_mat] = DataFile(LINKSDat,fibre1Dat,fibre2Dat,matrixDat);

%% Create LINKS datafile
LINKS_DatFile (RealizationFolder,MeshFolder,s_info,s_mat)

end

