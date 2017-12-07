function [lb, ub, IntCon, fibre1Dat, fibre2Dat, matrixDat, LINKSDat, ...
    RealizationFolder, MeshFolder, FolderName, s_info] = InitRVE(AdmissibleSet, ...
    PropStruct, GroupA, GroupB)


%% Bound constraints for the GA optimization problem
[sizeA,~]=size(PropStruct.Reference(PropStruct.Group==GroupA));
[sizeB,~]=size(PropStruct.Reference(PropStruct.Group==GroupB));
lb = [1, 1];           % Lower bound
ub = [sizeA, sizeB];   % Upper bound
IntCon = 1:2;          % Statement of the integer variables

%% Random initial fibre pair (Fn) generation
Fn(1) = randi([lb(1),ub(1)],1);
Fn(2) = randi([lb(2),ub(2)],1);

%% Mapping of the random chosen fibres to actual admissible reference numbers
[Fn] = MapVariables(Fn, AdmissibleSet);

%% Write fibreX.mat (X=1,2) - Fibre datafiles with fibre properties
WriteFibreDataFile(Fn,PropStruct);

%% RVE and mesh generation (chopped version of GenHybridComp.m)
[fibre1Dat,fibre2Dat,matrixDat,LINKSDat,RealizationFolder,MeshFolder,FolderName,s_info]=PREGenHybridComp();


end

