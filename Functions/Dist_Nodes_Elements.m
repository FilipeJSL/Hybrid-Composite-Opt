%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Function to distinguish node and element groups               %
%                                                                         %
%               Created by Lars Peperkamp October 2016 ©                  %
%                                                                         %
%              Based on function "f_mesh_quad_per_3D.m" by                %
%               Antonio Rui Melro - antonio.melro@fe.up.pt                %
%                             February 2009                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       10/2016 - First Version                                           %
%       10/2017 - Extended to more types of elements (Igor A.R.Lopes)     %
%       10/2017 - Particularised for voids and inclusions (Igor A.R.Lopes)%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nodes,Elements] = Dist_Nodes_Elements(Nodes,Elements,FolderPath)
disp(' ')
disp('Distinguishing Node and Element Groups')

%% Read Fibre positions
% Find file with Inclusions positions
FilePath = strcat(FolderPath,'\Fibre_Positions.txt');

% Read file
fileID = fopen(FilePath,'r');
formatSpec = '%f %f %f %f %f %f';
tline = 0;
while tline ~= -1 % Read lines until end of file is reached (tline = -1)
    FibrePos = textscan(fileID,formatSpec,'Delimiter','\n',...
        'CollectOutput', true);
    tline = fgetl(fileID);
end
fclose(fileID);

% Save data
FibrePos = FibrePos{1};
NumFibre = length(FibrePos);

%% Distinguish fibre nodes from matrix nodes
NodeMaterial (1:length(Nodes)) = 0;     % Vector for node material ID
for i = 1:length(Nodes)
    x1 = Nodes(i,2); y1 = Nodes(i,3);   % Read x and y position of nodes
    for j=1:NumFibre
        
        % Read x and y position of fibres
        xx = FibrePos(j,2); yy = FibrePos(j,3);
        DIST_TMP = sqrt((xx-x1)^2 + (yy-y1)^2);
        
        % If difference in x and y positions is small enough nodes are
        % part of fibre j
        if DIST_TMP <= 1.001*FibrePos(j,6)
            NodeMaterial(i) = j; % Node i belongs to fibre j
            break
        end
    end
end

%% Distinguish fibre elements from matrix elements
ngroups = 0;
if isfield(Elements,'ThreeNode')
    Elem3 = Elements.ThreeNode;
    Elem3(:,3:5) = Elem3(:,2:4);    % Add column to element matrix
    Elem3(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem3)           % Loop over elements
        k = 0;
        for j=3:5                   % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem3(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem3(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem3(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem3(i,j)) > 0 && NodeMaterial(Elem3(i,j)) == ...
                        NodeMaterial(Elem3(i,j-1))   % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 3               % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem3(i,2) = ft;    % Fibre type
        end
    end
    
    ngroups = max(Elem3(:,2));  % Save the number of elements groups found
    
    Elements.ThreeNode = Elem3; % Update Nodes Information
end

if isfield(Elements,'FourNode')
    Elem4 = Elements.FourNode;
    Elem4(:,3:6) = Elem4(:,2:5);    % Add column to element matrix
    Elem4(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem4)           % Loop over elements
        k = 0;
        for j=3:6                   % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem4(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem4(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem4(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem4(i,j)) > 0 && NodeMaterial(Elem4(i,j)) == ...
                        NodeMaterial(Elem4(i,j-1)) % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 4               % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem4(i,2) = ft;    % Fibre type
        end
    end
    
    Elem4(:,2) = Elem4(:,2) + ngroups; % Update the elements groups IDs
    ngroups = max(Elem4(:,2));         % and save their number
    
    Elements.FourNode = Elem4; % Update Nodes Information
end

if isfield(Elements,'SixNode')
    Elem6 = Elements.SixNode;
    Elem6(:,3:8) = Elem6(:,2:7);    % Add column to element matrix
    Elem6(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem6)           % Loop over elements
        k = 0;
        for j=3:8                   % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem6(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem6(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem6(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem6(i,j)) > 0 && NodeMaterial(Elem6(i,j)) == ...
                        NodeMaterial(Elem6(i,j-1)) % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 6               % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem6(i,2) = ft;    % Fibre type
        end
    end
    
    Elem6(:,2) = Elem6(:,2) + ngroups; % Update the elements groups IDs
    ngroups = max(Elem6(:,2));         % and save their number
    
    Elements.SixNode = Elem6; % Update Nodes Information
end

if isfield(Elements,'EightNode')
    Elem8 = Elements.EightNode;
    Elem8(:,3:10) = Elem8(:,2:9);   % Add column to element matrix
    Elem8(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem8)           % Loop over elements
        k = 0;
        for j=3:10                  % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem8(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem8(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem8(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem8(i,j)) > 0 && NodeMaterial(Elem8(i,j)) == ...
                        NodeMaterial(Elem8(i,j-1)) % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 8               % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem8(i,2) = ft;    % Fibre type
        end
    end
    
    Elem8(:,2) = Elem8(:,2) + ngroups; % Update the elements groups IDs
    ngroups = max(Elem8(:,2));         % and save their number
    
    Elements.EightNode = Elem8; % Update Nodes Information
end

if isfield(Elements,'TenNode')
    Elem10 = Elements.TenNode;
    Elem10(:,3:12) = Elem10(:,2:11);% Add column to element matrix
    Elem10(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem10)          % Loop over elements
        k = 0;
        for j=3:12                  % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem10(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem10(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem10(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem10(i,j)) > 0 && NodeMaterial(Elem10(i,j)) == ...
                        NodeMaterial(Elem10(i,j-1)) % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 10              % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem10(i,2) = ft;   % Fibre type
        end
    end
    
    Elem10(:,2) = Elem10(:,2) + ngroups; % Update the elements groups IDs
    ngroups = max(Elem10(:,2));          % and save their number
    
    Elements.TenNode = Elem10; % Update Nodes Information
end

if isfield(Elements,'TwentyNode')
    Elem20 = Elements.TwentyNode;
    Elem20(:,3:22) = Elem20(:,2:21);% Add column to element matrix
    Elem20(:,2) = 3;                 % Change value of new column to 3, because
                                    % Fibre1 = 1, Fibre2 = 2, Matrix = 3
    for i=1:length(Elem20)          % Loop over elements
        k = 0;
        for j=3:22                  % Loop over element nodes
            if j == 3
                if NodeMaterial(Elem20(i,j)) > 0      % If node is a fibre node
                    k = k + 1;
                    if FibrePos(NodeMaterial(Elem20(i,j)),5)==0   % Fibre1 Node
                        ft = 1;
                    elseif FibrePos(NodeMaterial(Elem20(i,j)),5)==1   % Fibre2 Node
                        ft = 2;
                    end
                end
            else
                if NodeMaterial(Elem20(i,j)) > 0 && NodeMaterial(Elem20(i,j)) == ...
                        NodeMaterial(Elem20(i,j-1)) % Check other element nodes
                    k = k + 1;
                end
            end
        end
        if k == 20              % If all elements nodes are fibre nodes
                                % add fibre type to element
            Elem20(i,2) = ft;   % Fibre type
        end
    end
    
    Elem20(:,2) = Elem20(:,2) + ngroups; % Update the elements groups IDs
    
    Elements.TwentyNode = Elem20; % Update Nodes Information
end
end