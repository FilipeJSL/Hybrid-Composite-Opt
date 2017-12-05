function [] = PREGenHybridComp()

%% Global variables
global Abaqus elem_size_param

%% Add folders to Matlab path and define input files
addpath('.\Functions'); % Folder with all functions
RVEDat = 'Datafiles\RVE.mat ';       % File with the composite RVE properties
fibre1Dat = 'Datafiles\fibre1.mat '; % File with the fibre type 1 properties
fibre2Dat = 'Datafiles\fibre2.mat '; % File with the fibre type 2 properties
LINKSDat  = 'Datafiles\LINKSOptions.mat'; % File with LINKS datafile parameters

Abaqus    = 'C:\SIMULIA\Commands\abaqus'; % Abaqus executable

%% Define the element size parameter
elem_size_param = 2.5;

%% Read  RVE inputs
[Input] = RVEinput(RVEDat,fibre1Dat,fibre2Dat);

%% Read LINKS input options
[s_info,~] = DataFile(LINKSDat);

%% RVE Generation
% Name of output folder
FolderName = 'Output';
RealizationFolder = strcat(FolderName,'\realization_','StandardRVE');

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

end

