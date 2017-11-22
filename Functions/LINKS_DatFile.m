%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Function to generate LINKS input file                  %
%                                                                         %
%             Created by Igor A. Rodrigues Lopes, October 2017            %
%              Based on the function Input.m (Lars Peperkamp)             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       10/2017 - First Version                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Used Functions:                                                       %
%                                                                         %
%       Nodes_Elements.m        - Function to read nodes and elements     %
%                                                                         %
%       Dist_Nodes_Elements.m   - Function to distinguish node and        %
%                                 element groups                          %
%       DataFile.m              - Function to read data: title,           %
%                                 analysis type, materials etc.           %
%       Write_Output.m          - Function to write the output file       %
%                                                                         %
%       SubSubFolder.m          - Function to find paths to subsubfolders %  
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LINKS_DatFile(RealizationFolder, MeshFolder, s_info,s_mat)
%%
disp(' ')
disp('Reading Mesh File')

%% Read nodes and elements
MeshInp = strcat(MeshFolder,'\mesh.inp');
fileID = fopen(MeshInp,'r');  % Open Abaqus mesh file

if s_info.ntype == 6
    NDim = 3;
else
    NDim = 2;
end

[Nodes,Elements] = Nodes_Elements(fileID,NDim);

fclose(fileID);                 % Close Abaqus mesh file

%% Distinquish Matrix and Fibre elements
[Nodes,Elements] = Dist_Nodes_Elements(Nodes,Elements,RealizationFolder);

%% Save output file .dat
disp(' ')
disp('Writing Data')

% There are multiple Boundary options?
for ibound = 1:length(s_info.BC)
    % Create datafile for each boundary condition
    s_info.bound_condition = strcat(s_info.BC(ibound),'_Condition');
    Full_Path = strcat(MeshFolder,'\RVE_',s_info.BC(ibound),'.rve');
    Write_Output(Full_Path,s_info,s_mat,Nodes,Elements);
    Full_Path = {Full_Path}; % Change path to string to make it compatible
                             % with multiple input files
end

disp(' ')
disp('Finished - Writing Data')
disp(' ')
disp('LINKS Input file(s) stored in')
disp(MeshFolder)

end