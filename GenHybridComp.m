%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         Script to generate a random RVE for hybrid composites           %
%                                                                         %
%                        Igor A. Rodrigues Lopes                          %
%                            ilopes@fe.up.pt                              %
%                                                                         %
%                             November 2017                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       11/2017 - Initial coding. Based on MSP_Random_RVE                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This is coded to generate a RVE with randomly distributed circular      %
% inclusions and voids. The input is:                                     %
%   - the total volume fraction of inclusions and voids                   %
%   - the total number of inclusions and voids                            %
%   - the number of voids                                                 %
%   - the radius of the inclusions and voids (it is the same)             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear variables
clear all
close all
clc

%% Global variables
global Abaqus elem_size_param
%% Add folders to Matlab path and define input files
addpath('.\Functions'); % Folder with all functions
RVEDat = 'Datafiles\RVE.mat ';       % File with the composite RVE properties
fibre1Dat = 'Datafiles\fibre1.mat '; % File with the fibre type 1 properties
fibre2Dat = 'Datafiles\fibre2.mat '; % File with the fibre type 2 properties
matrixDat = 'Datafiles\matrix.mat '; % File with the matrix properties
LINKSDat  = 'Datafiles\LINKSOptions.mat'; % File with LINKS datafile parameters
Abaqus    = 'C:\SIMULIA\Commands\abaqus'; % Abaqus executable

%% Define the number of realizations to be generated
nrlz = 1;    % number of realizations

%% Define the element size parameter
elem_size_param = 3.0;

%% Read  RVE inputs
[Input] = RVEinput(RVEDat,fibre1Dat,fibre2Dat);

%% Read LINKS input options
[s_info,s_mat] = DataFile(LINKSDat);

%% Loop over realizations
for irlz = 1:nrlz
    %% Name of output folder
    FolderName = 'Output';
    RealizationFolder = strcat(FolderName,'\realization_',num2str(irlz));

    mkdir(RealizationFolder)

    %% Generate Fibre Distribution
    Create_Random_RVE (Input, RealizationFolder);

    %% Build Mesh (with Abaqus)
    MeshFolder = strcat(RealizationFolder,'\mesh');
    mkdir(MeshFolder)

    % Define dimension
    if s_info.ntype == 6
        NDim = 3;
    else
        NDim = 2;
    end

    status_mesh = 0;
    while status_mesh == 0
        % Run function to generate mesh
        if(NDim == 2)
            status_mesh = Create_2D_Mesh(MeshFolder);
        elseif(NDim == 3)
            status_mesh = Create_3D_Mesh(MeshFolder);
        else
            error('Invalid dimension')
        end
    end

    %% And create LINKS datafile
    LINKS_DatFile (RealizationFolder, MeshFolder, s_info,s_mat)
end

%% Run LINKS for the generated RVEs

% add here some commands to create a run_LINKS.bat

%% Plot the stress evolution and 

% Define the path for the .out file of the realization and call the
% plotstress function

% for irlz = 1:nrlz
%     plotstress(LINKSout_path, irlz);
% end