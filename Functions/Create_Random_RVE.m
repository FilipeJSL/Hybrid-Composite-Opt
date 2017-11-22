%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%            Function that creates a random RVE for hybrid composites     %
%                                                                         %
%                 Igor A. Rodrigues Lopes, November 2017                  %
%                              ilopes@fe.up.pt                            %
%                                                                         %
%            Based on Matlab routine "Create_Random_RVE_Mesh.m"           %
%                created by Lars Peperkamp October 2016                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       11/2017 - Initial coding                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Used Functions:                                                       %
%                                                                         %
%         RAND_PER_uSTRU_GEN_3D - Generates a random distribution of      %
%                                 circles that represent the real fibre   %
%                                 distribution                            %
%         f_image_per           - Generates images of fibre distribution  %
%                                 in .BMP                                 %
%         Create_2D_Mesh        - Function that writes Python script to   %
%                                 generate a 2D mesh in Abaqus            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Create_Random_RVE(Input, OutputFolder)
tic; % Set counter

%% Definition of Global Variables
global R delta deltaoop Vol_fibre_req DISTMIN Max_fibres a b c ...
    N_guesses_max N_cycles_max N_change Square_size Square_inc ...
    S_base Fibre_type_1 R1 Rmax Rmin Fibre_pos elem_size_param;

%% Display message
disp(' '); disp('Setting Random RVE Parameters')

%% Set Fibres Diameters
R  = Input.Radius1;   % Fibre1 radius
R1 = Input.Radius2;   % Fibre2 radius

%%  Set Total Volume Fraction

Vol_fibre_req = Input.TotalVolFrac; % Total fibre volume

%% Set relative fibre volume fractions
Fibre_type_2 = Input.Fibre2VolFrac;
Fibre_type_1 = 1 - Fibre_type_2;

%% Set RVE size parameter
delta    = Input.RVEsize;  % RVE in-plane size
deltaoop = Input.RVEthick; % RVE out-of-plane size
%% Store Input Parameters
Rmax            =   max(R,R1);               % Maximal fibre radius
Rmin            =   min(R,R1);               % Minimal fibre radius
% Minimum distance between fibre centres %%%%% 2.07
DISTMIN         =   2.07;
% Maximum number of fibres
Max_fibres      =   750;
% Number of guesses in the random generation of fibre position
N_guesses_max   =   50E3;
% Maximum number of cycles that the routine runs
N_cycles_max    =   30;
% Number of iterations before changing criteria on First Heuristic
% MinValue = 3
N_change        =   3;
% Initial size of square in Second Heuristic
Square_size     =   3*Rmax;
% Increment to be given to the square size in Second Heuristic
Square_inc      =   (8.5-10*Vol_fibre_req)*Rmax;
% Average of element size
S_base          =   1.5*2*pi/60*Rmin;

% I will double the element size since quadratic elements are used now (Igor)
S_base = elem_size_param*S_base;

%% Determination of RVE Dimensions
a = delta*(Rmax+Rmin)/2; % RVE width (yy direction)
b = delta*(Rmax+Rmin)/2; % RVE height (zz direction)
c = deltaoop*(Rmax+Rmin)/2; % RVE thickness (xx direction)

%% Generate Fibre Distribution
status = 0;
while status == 0
    status = RAND_PER_uSTRU_GEN_3D;
    if status == 1, break; end
end

%% Print image of fibre distribution
index1 = 0;
f_image_per(index1,OutputFolder);

%% Save fibres positions in .txt file
Print_txt(OutputFolder,'\Fibre_Positions.txt',Fibre_pos, ...
                            '%d %d %d %d %d %d \r\n');

end