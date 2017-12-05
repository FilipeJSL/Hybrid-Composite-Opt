function [fvalue] = OptObjectiveFun(Fn, AdmissibleSet, PropStruct)

%% Mapping of the GA design variables to the admissible set of fibres
Fn=MapVariables(Fn, AdmissibleSet);

%% Write fibreX.mat (X=1,2) - Fibre datafiles with fibre properties
WriteFibreDataFile(Fn,PropStruct);

%% Call to GenHybridComp.m
GenHybridComp

%%
end

