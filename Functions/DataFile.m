function [s_info,s_mat]=DataFile(LINKSDat)
%  U.PORTO-FEUP-DEMec
%  Thesis
%  Miguel Viera de Carvalho (em11129@fe.up.pt)(miguel.carvalho.129@gmail.com)
%  created: 20/06/2016
%  revised: 10/2016 by Lars Peperkamp
%           11/2016 by Lars Peperkamp for use with Von Mises model and
%           adaptable Deformation Gradients
%  modified: 10/2017 by Igor A. Rodrigues Lopes
%            Make it compatible with LINKS

%% Function to read data: title, analysis tipe, materials
% title, analysis type, deformation gradient, step length,
% max number of line searchs, number of increments,
% non local formulation, boundary coundition, material info


%% Store data from Main input file in cell array
disp(' '); disp('Reading LINKS Options')

formatSpec = '%s'; %Read file as string
fileID = fopen(LINKSDat,'r');

%Read file and store every line on own line in cell
Input = textscan(fileID,formatSpec,'Delimiter','\n');
fclose(fileID);

%% Title
% Combine main and material input to title for MSP input file
TitleStr = strfind(Input{1},'Title');
for i = 1:length(TitleStr)
    empty = isempty(cell2mat(TitleStr(i)));
    if empty == 0
        TitleLoc = i+1;
    end
end

title = Input{1}(TitleLoc);
title = title{1};

%% Analysis Type
AnaTypeStr = strfind(Input{1},'Analysis Type');
for i = 1:length(AnaTypeStr)
    empty = isempty(cell2mat(AnaTypeStr(i)));
    if empty == 0
        AnaTypeLoc = i+1;
    end
end

ntype = str2double(Input{1}(AnaTypeLoc));

%% Step Length Value
USLStr = strfind(Input{1},'Upper Step Length Value');
for i = 1:length(USLStr)
    empty = isempty(cell2mat(USLStr(i)));
    if empty == 0
        USLLoc = i+1;
    end
end

LSLStr = strfind(Input{1},'Lower Step Length Value');
for i = 1:length(LSLStr)
    empty = isempty(cell2mat(LSLStr(i)));
    if empty == 0
        LSLLoc = i+1;
    end
end

USL = str2double(Input{1}(USLLoc));
LSL = str2double(Input{1}(LSLLoc));

% column 1 -> lower
% column 2 -> upper
step_len = [LSL USL];

%% Maximum Number of Line Searchs
NLSStr = strfind(Input{1},'Maximum Number of Line Searches');
for i = 1:length(NLSStr)
    empty = isempty(cell2mat(NLSStr(i)));
    if empty == 0
        NLSLoc = i+1;
    end
end

max_search = str2double(Input{1}(NLSLoc));

%% Number of Increments
NOIStr = strfind(Input{1},'Number of Increments');
for i = 1:length(NOIStr)
    empty = isempty(cell2mat(NOIStr(i)));
    if empty == 0
        NOILoc = i+1;
    end
end

num_incr = str2double(Input{1}(NOILoc));

%% Convergence tolerance
NOIStr = strfind(Input{1},'Convergence tolerance');
for i = 1:length(NOIStr)
    empty = isempty(cell2mat(NOIStr(i)));
    if empty == 0
        NOILoc = i+1;
    end
end

conv_tol = str2double(Input{1}(NOILoc));

%% Maximum Number of Increments Cuts
NOIStr = strfind(Input{1},'Maximum number of consecutive increment cuts');
for i = 1:length(NOIStr)
    empty = isempty(cell2mat(NOIStr(i)));
    if empty == 0
        NOILoc = i+1;
    end
end

max_cut = str2double(Input{1}(NOILoc));

%% Solvers
SolStr = strfind(Input{1},'* Solver');
for i = 1:length(SolStr)
    empty = isempty(cell2mat(SolStr(i)));
    if empty == 0
        SolLoc = i+1;
    end
end

SolParStr = strfind(Input{1},'Parallel Solver');
for i = 1:length(SolParStr)
    empty = isempty(cell2mat(SolParStr(i)));
    if empty == 0
        SolParLoc = i+1;
    end
end

solver = Input{1}(SolLoc);
parallel_solver = str2double(Input{1}(SolParLoc));

%% VTK output
VTKStr = strfind(Input{1},'VTK output');
for i = 1:length(VTKStr)
    empty = isempty(cell2mat(VTKStr(i)));
    if empty == 0
        VTKLoc = i+1;
    end
end

VTKout = Input{1}(VTKLoc);

%% Arc Length
% options:  % ON
% OFF
arc_length = 'OFF';         %Don't change feature is not well established

%% Type of coordinate system
% options:  % CARTESIAN
% CYLINDRICAL
coord_system = 'CARTESIAN'; % Don't change MSP only works with Cartesian

%% Non Local Formulation
% options:  % ON
% OFF
n_local_for = 'OFF';        % Don't change MSP_Random_RVE needs to be
% adapted to use this function

%% Assign to struct variable
s_info.title=title;
s_info.coord_system=coord_system;
s_info.ntype = ntype;
s_info.arc_length = arc_length;
s_info.step_len = step_len;
s_info.max_search = max_search;
s_info.num_incr = num_incr;
s_info.conv_tol = conv_tol;
s_info.max_cut = max_cut;
s_info.n_local_for = n_local_for;
s_info.solver = solver{1};
s_info.parallel_solver = parallel_solver;
s_info.VTKout = VTKout;

%% Stretch Ratio
StretchStr = strfind(Input{1},'Stretch Ratio');
for i = 1:length(StretchStr)
    empty = isempty(cell2mat(StretchStr(i)));
    if empty == 0
        StretchLoc = i+1;
    end
end

Stretch = str2double(Input{1}(StretchLoc));

s_info.Stretch = Stretch;


%% Read Boundary Conditions
BCStr = strfind(Input{1},'* Boundary Condition');
for i = 1:length(BCStr)
    empty = isempty(cell2mat(BCStr(i)));
    if empty == 0
        BCLoc = i+1;
    end
end

bound_condition = Input{1}(BCLoc);

% Split the strings and remove any empty strings.
BC = strsplit(bound_condition{1},'_Condition');
LBC = length(BC);
for i = 1:LBC
    if isempty(BC{i}) == 1
        if i == 1
            BC = BC(2:length(BC));
            break
        elseif i == LBC
            BC = BC(1:(length(BC)-1));
            break
        end
    end
end
s_info.BC = strtrim(BC);


%% Check Multiple F or BC Options   
FStr = strfind(Input{1},'Prescribed Deformation Gradient');
for i = 1:length(FStr)
    empty = isempty(cell2mat(FStr(i)));
    if empty == 0
        FLoc = i+1;
    end
end
    
F = str2double(Input{1}(FLoc));

if ntype == 6
    ndim = 3;
else
    ndim = 2;
end

def_grad = eye(ndim);
    
if F == 1
    def_grad(1,1) = Stretch;
elseif F == 2
    def_grad(2,2) = Stretch;
elseif F == 3
    def_grad(3,3) = Stretch;
elseif F == 12
    def_grad(1,1) = Stretch-1;
elseif F == 112
    def_grad(1,1) = Stretch;
    def_grad(1,2) = Stretch-1;
elseif F == 212
    def_grad(2,2) = Stretch;
    def_grad(1,2) = Stretch-1;
elseif F == 21
    def_grad(1,1) = Stretch;
    def_grad(2,2) = Stretch;
elseif F == 2112
    def_grad(2,1) = Stretch-1;
    def_grad(1,2) = Stretch-1;
end

s_info.def_grad = def_grad;
s_info.F = F;

%% Read file with material properties
disp(' ')
disp('Reading Material Properties')

% Three distinct material types: 1 - fibre 1, 2 - fibre 2, 3 - matrix
s_mat.nmats = 3;
materials = {'Datafiles\fibre1.mat' 'Datafiles\fibre2.mat' ...
             'Datafiles\matrix.mat'};

formatSpec = '%s'; %Read file as string
for imaterial = 1:s_mat.nmats
    % Open file with material properties
    fileID_mat = fopen(cell2mat(materials(imaterial)),'r');
    % Read file and put every line on own line in cell
    A = textscan(fileID_mat,formatSpec,'Delimiter','\n');
    fclose(fileID_mat);

    % Material type
    MatTypStr = strfind(A{1},'Material type');
    for i = 1:length(MatTypStr)
        empty = isempty(cell2mat(MatTypStr(i)));
        if empty == 0
            MatTypLoc = i+1;
        end
    end

    mat_type = A{1}(MatTypLoc);

    % Material density
    MatTypStr = strfind(A{1},'Density');
    for i = 1:length(MatTypStr)
        empty = isempty(cell2mat(MatTypStr(i)));
        if empty == 0
            MatTypLoc = i+1;
        end
    end

    mat_dense = A{1}(MatTypLoc);

    % Material properties (inline)
    MatTypStr = strfind(A{1},'List of properties');
    for i = 1:length(MatTypStr)
        empty = isempty(cell2mat(MatTypStr(i)));
        if empty == 0
            MatTypLoc = i+1;
        end
    end

    mat_props = A{1}(MatTypLoc);

    % Number of hardening points
    MatTypStr = strfind(A{1},'Number of Hardening points');
    for i = 1:length(MatTypStr)
        empty = isempty(cell2mat(MatTypStr(i)));
        if empty == 0
            MatTypLoc = i+1;
        end
    end

    mat_nhard = str2double(A{1}(MatTypLoc));

    % List of hardening points
    MatTypStr = strfind(A{1},'List of Hardening points');
    for i = 1:length(MatTypStr)
        empty = isempty(cell2mat(MatTypStr(i)));
        if empty == 0
            MatTypLoc = i+1;
        end
    end

    for i = 1:mat_nhard
        mat_hardp(i) = A{1}(MatTypLoc);
        MatTypLoc = MatTypLoc + 1;
    end

    % store in structured variable
    s_mat.type(imaterial)  = mat_type;
    s_mat.dense(imaterial) = mat_dense;
    s_mat.props(imaterial) = mat_props;
    s_mat.nhard(imaterial) = mat_nhard;
    if mat_nhard ~= 0
        s_mat.hardp(:,imaterial) = mat_hardp;
    end
end

clearvars -except s_info s_mat

end