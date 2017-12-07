function [fvalue] = OptObjectiveFun(Fn, AdmissibleSet, PropStruct, fibre1Dat, ...
    fibre2Dat, matrixDat, LINKSDat, RealizationFolder, MeshFolder, FolderName, ...
    CurrentFolder, s_info)

%% Mapping of the GA design variables to the admissible set of fibres
Fn=MapVariables(Fn, AdmissibleSet);

%% Write fibreX.mat (X=1,2) - Fibre datafiles with fibre properties
WriteFibreDataFile(Fn,PropStruct);

%% LINKS input datafile generation (chopped version of GenHybridComp.m)
POSTGenHybridComp (LINKSDat,fibre1Dat,fibre2Dat,matrixDat,RealizationFolder,MeshFolder,s_info);

%% Run LINKS for the standard RVE with the set of fibres defined by Fn
[rveFile]=RunLINKS(MeshFolder, FolderName, CurrentFolder, s_info);

%% Compute the objective value
[E] = ComputeObj(MeshFolder, rveFile, s_info);

fprintf('Fibre 1: %f  Fibre 2: %f  Stiffness: %g',Fn(1),Fn(2),E);

fvalue=-E;


end

