%% Clear variables
clear all
close all
clc

tic

%% Add optimization folder to Matlab path
addpath('.\Functions Opt'); % Folder with all optimization related functions
FibreMatDat  = 'Datafiles\FibreMaterialData.mat'; % Database file with all fibres' properties

%% Load fibre properties database file
[PropStruct] = GetFibreProps(FibreMatDat);

%% Fibre groups to composite optimization (see FibreMaterialData.mat)
GroupA=1;
GroupB=3;

AdmSet1=PropStruct.Reference(PropStruct.Group==GroupA);
AdmSet2=PropStruct.Reference(PropStruct.Group==GroupB);
AdmissibleSet={AdmSet1 AdmSet2};

%% RVE Initialization (This module keeps the RVE fixed throughout the optimization process)
[lb,ub,IntCon] = InitRVE(AdmissibleSet,PropStruct,GroupA,GroupB);

%% MATLAB's GA tool setup
fitnessFunction = @(Fn)OptObjectiveFun(Fn, AdmissibleSet, PropStruct);  % Function handle to the fitness function
numberOfVariables = 2;   % Number of decision variables

%Options to refine/loose the search for the optimal solution (there are other parameters that can be set)
PopulationSize = 150;
MaxGenerations = 200;
EliteCount = 10;
FunctionTolerance = 1e-8;

options = optimoptions(@ga, ...
                    'PopulationSize', PopulationSize, ...
                    'MaxGenerations', MaxGenerations, ...
                    'EliteCount', EliteCount, ...
                    'FunctionTolerance', FunctionTolerance, ...
                    'PlotFcn', @gaplotbestf);

% Call to the GA tool
[xbest, fbest, exitflag] = ga(fitnessFunction, numberOfVariables, [], [],...
    [], [], lb, ub, [], IntCon, options);

%%
ElapsedTime=toc;
fprintf('Total optimization time: %.2f minutes\n',ElapsedTime/60);
