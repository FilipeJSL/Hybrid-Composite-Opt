function [] = Write_Output(sv_Full_Path,s_info,s_mat,Nodes,Elements)
%  U.PORTO-FEUP-DEMec
%  Thesis
%  Miguel Viera de Carvalho (em11129@fe.up.pt)(miguel.carvalho.129@gmail.com)
%  created: 29/06/2016
%  revised: 10/2016 by Lars Peperkamp for use in MSP_Random_RVE
%           11/2016 by Lars Peperkamp for use with Von Mises model
%%
global c
%% Function to write the output file

%% DataFile Info
title = s_info.title;
ntype = s_info.ntype;
coord_system = s_info.coord_system;
def_grad = s_info.def_grad;
step_len = s_info.step_len;
max_search = s_info.max_search;
num_incr = s_info.num_incr;
max_cut = s_info.max_cut;
conv_tol = s_info.conv_tol;
n_local_for = s_info.n_local_for;
bound_condition = s_info.bound_condition;
solver = s_info.solver;
parallel_solver = s_info.parallel_solver;
VTKout = s_info.VTKout{1};

%% Element Groups
SizeEleGroup = 0;
SizeEleType  = 0;
SizeEle = 0;

if isfield(Elements,'ThreeNode')
    Ele3 = Elements.ThreeNode;
    SizeEleType = SizeEleType + 1;
    EleName{SizeEleType} = 'TRI_3';
    NGP(SizeEleType) = 1;
    for i = SizeEleGroup + 1 : max(Ele3(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele3(:,2));
    SizeEle = SizeEle + length(Ele3);
end

if isfield(Elements,'FourNode')
    Ele4 = Elements.FourNode;
    SizeEleType = SizeEleType + 1;
    if ntype == 6
        EleName{SizeEleType} = 'TETRA_4';
        NGP(SizeEleType) = 1;
    else
        EleName{SizeEleType} = 'QUAD_4';
        NGP(SizeEleType) = 4;
    end
    for i = SizeEleGroup + 1 : max(Ele4(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele4(:,2));
    SizeEle = SizeEle + length(Ele4);
end

if isfield(Elements,'SixNode')
    Ele6 = Elements.SixNode;
    SizeEleType = SizeEleType + 1;
    EleName{SizeEleType} = 'TRI_6';
    NGP(SizeEleType) = 1;
    for i = SizeEleGroup + 1 : max(Ele6(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele6(:,2));
    SizeEle = SizeEle + length(Ele6);
end

if isfield(Elements,'EightNode')
    Ele8 = Elements.EightNode;
    SizeEleType = SizeEleType + 1;
    if ntype == 6
        EleName{SizeEleType} = 'HEXA_8';
        NGP(SizeEleType) = 8;
    else
        EleName{SizeEleType} = 'QUAD_8';
        NGP(SizeEleType) = 4;
    end
    for i = SizeEleGroup + 1 : max(Ele8(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele8(:,2));
    SizeEle = SizeEle + length(Ele8);
end

if isfield(Elements,'TenNode')
    Ele10 = Elements.TenNode;
    SizeEleType = SizeEleType + 1;
    EleName{SizeEleType} = 'TETRA_10';
    NGP(SizeEleType) = 4;
    for i = SizeEleGroup + 1 : max(Ele10(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele10(:,2));
    SizeEle = SizeEle + length(Ele10);
end

if isfield(Elements,'TwentyNode')
    Ele20 = Elements.TwentyNode;
    SizeEleType = SizeEleType + 1;
    EleName{SizeEleType} = 'HEXA_20';
    NGP(SizeEleType) = 8;
    for i = SizeEleGroup + 1 : max(Ele20(:,2))
        Ele_Type(i) = SizeEleType;
        Ele_Mat(i) = i - SizeEleGroup;
    end
    SizeEleGroup = max(Ele20(:,2));
    SizeEle = SizeEle + length(Ele20);
end

%% Write output file
% Output File ID
sv_Full_Path = sv_Full_Path{1};
ID = fopen(sv_Full_Path,'w');

fprintf(ID,'%s\r\n%s\r\n', 'TITLE', title);
fprintf(ID,'\r\n');
fprintf(ID,'%s\r\n', 'Prescribed_Deformation_Gradient'); 
if ntype == 6
    formatF = '%g %g %g \r\n';
else
    formatF = '%g %g \r\n';
end
fprintf(ID,formatF, def_grad');
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'ANALYSIS_TYPE', ntype);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %s\r\n', 'SOLVER',solver);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'PARALLEL_SOLVER', parallel_solver);
fprintf(ID,'\r\n');
fprintf(ID,'%s\r\n%g\r\n', 'STEP_LENGTH_LOWER_VALUE',step_len(1));
fprintf(ID,'\r\n');
fprintf(ID,'%s\r\n%g\r\n', 'STEP_LENGTH_UPPER_VALUE',step_len(2));
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'MAXIMUM_NUMBER_LINE_SEARCHS',max_search);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'Number_of_Increments', num_incr);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'MAX_CONSECUTIVE_INCREMENT_CUTS', max_cut);
fprintf(ID,'\r\n');
fprintf(ID,'%s\r\n%g\r\n', 'CONVERGENCE_TOLERANCE', conv_tol);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %s\r\n', 'NON_LOCAL_FORMULATION', n_local_for);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %s\r\n', 'VTK_OUTPUT', VTKout);
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'ELEMENT_GROUPS',SizeEleGroup); % number of groups
for i = 1:SizeEleGroup
    % write the matrix with:
    % column 1: group number
    % column 2: element type number
    % column 3: material number
    fprintf(ID,'    %g  %g    %g\r\n', i, Ele_Type(i),Ele_Mat(i));
end
fprintf(ID,'\r\n');
fprintf(ID,'%s  %g\r\n', 'ELEMENT_TYPES', SizeEleType); % number of types
for i = 1:SizeEleType
    fprintf(ID,'    %g  %s\r\n ', i,EleName{i});
    if strcmp(EleName{i},'TRI_3') == 0       % TRI_3 does not have the Gauss
        fprintf(ID,'        %g  %s\r\n',...  % points number defined
            NGP(i),'GP');
    end
end
fprintf(ID,'\r\n');
SizeMat = s_mat(1).nmats;
fprintf(ID,'%s  %g\r\n', 'MATERIALS', SizeMat);
for i = 1:SizeMat % write material info
    type = s_mat.type(i);
    dense = s_mat.dense(i);
    if (type{1}(1:4) == 'UMAT')
        fprintf(ID,'    %g  %s\r\n', i, type{1});
        % write the RVE length to the last property if is UMATFIBRAND
        if (type{1}(1:11) == 'UMATFIBRAND')
            s_mat.hardp{12,i} = [num2str(c),' ',s_mat.hardp{12,i}];
        end
    else
        props = s_mat.props(i);
        fprintf(ID,'    %g  %s\r\n  %s \r\n  %s\r\n', i, type{1},...
            dense{1},props{1});
    end
    if s_mat.nhard(i) ~= 0
        fprintf(ID,'  %g\r\n',s_mat.nhard(i));
        for j = 1:s_mat.nhard(i)
            fprintf(ID,'  %s \r\n', s_mat.hardp{j,i}');
        end
    end
end
fprintf(ID,'\r\n');
fprintf(ID,'%s  %s\r\n', 'Boundary_Type', bound_condition{1});
fprintf(ID,'\r\n');
fprintf(ID,'%s  %d  %s\r\n', 'NODE_COORDINATES', size(Nodes,1), coord_system);
if ntype == 6 % 3D
    % change coordinates sort, for fibres in x direction
    fprintf(ID,'    %d  %e  %e %e \r\n', Nodes(:,[1 4 2 3])');
else          % 2D
    fprintf(ID,'    %d  %e  %e \r\n', Nodes');
end
fprintf(ID,'\r\n');
fprintf(ID,'%s  %d \r\n', 'ELEMENTS', SizeEle);
fprintf(ID,'\r\n');
if isfield(Elements,'ThreeNode')
    for i = 1:length(Ele3)
        fprintf(ID,'    %d  %d  %d  %d  %d \r\n',Ele3(i,:)');
    end
end
if isfield(Elements,'FourNode')
    for i = 1:length(Ele4)
        fprintf(ID,'    %d  %d  %d  %d  %d  %d \r\n',Ele4(i,:)');
    end
end
if isfield(Elements,'SixNode')
    for i = 1:length(Ele6)
        fprintf(ID,'    %d  %d  %d  %d  %d  %d  %d  %d \r\n',Ele6(i,:)');
    end
end
if isfield(Elements,'EightNode')
    for i = 1:length(Ele8)
        fprintf(ID,'    %d  %d  %d  %d  %d  %d %d  %d  %d  %d \r\n',Ele8(i,:)');
    end
end
if isfield(Elements,'TenNode')
    for i = 1:length(Ele10)
        fprintf(ID,'    %d  %d  %d  %d  %d  %d %d  %d  %d  %d  %d  %d \r\n',Ele10(i,:)');
    end
end
if isfield(Elements,'TwentyNode')
    for i = 1:length(Ele20)
        fprintf(ID,'    %d  %d  %d  %d  %d  %d %d  %d  %d  %d  %d  %d  %d  %d  %d  %d %d  %d  %d  %d  %d  %d \r\n',Ele20(i,:)');
    end
end
fclose(ID);
end