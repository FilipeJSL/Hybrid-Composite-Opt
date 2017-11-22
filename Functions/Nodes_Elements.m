%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Function to read nodes and elements                    %
%                                                                         %
%                Created by Lars Peperkamp October 2016 ©                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       10/2016 - First Version                                           %
%       10/2017 - Extended to any element type (Igor A.R. Lopes)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nodes,Elements]=Nodes_Elements(fileID,NDim)
%% Read Nodes
disp(' '); disp('Reading Nodes')

filepos = 0;
tline = fgetl(fileID);

while tline ~= -1 % Read lines until end of file is reached (tline = -1)
    
    % When the line starts with *Node break while loop and save position
    % in file.
    if strncmp(tline,'*Node',5)
        filepos = ftell(fileID); % Position in file
        break % End while loop
    end
    
    filepos = ftell(fileID);
    tline = fgetl(fileID);
end

if tline == -1
    disp(' ')
    error('End of file is reached. No nodes and appropriate coordinates found.')
else
    % Go to the beginning of first line with nodes
    fseek(fileID, filepos, 'bof');
    
    % Use textscan to read nodes
    if(NDim == 2)
        formatSpec = '%f %f %f';
    elseif(NDim == 3)
        formatSpec = '%f %f %f %f';
    end
    Nodes = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
    filepos = ftell(fileID);
    
    % Convert cell array to array of doubles or integers
    Nodes = Nodes{1};
end

%% Read Elements
disp(' '); disp('Reading Elements')

nelemtyp = 0;
while nelemtyp ~= -1
% Read first type of elements
    while tline ~= -1 % Read lines until end of file is reached (tline = -1)
        % When the line starts with *Element break while loop and save position
        if strncmp(tline,'*Element',8);
            nelemtyp = nelemtyp + 1;
            filepos = ftell(fileID); % Position in file
            break % End tline while loop
        end
        filepos = ftell(fileID);
        tline = fgetl(fileID);
    end

    if tline == -1
        if nelemtyp == 0
            error('End of file is reached. No elements found.')
        else
            disp('End of file is reached. All element types have been read.')
            break % End nelemtyp while loop
        end
    else
        if isempty(strfind(tline,'CPS3')) == 0 % Check elemenent type
            if NDim ~= 2
                error('CPS3 element can only be used in 2D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d'; % Format of 3-Node elements: 4 doubles
            Ele3Node = textscan(fileID,formatSpec,'Delimiter',',',...
                'CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele3Node = Ele3Node{1};
            % Store in structure
            Elements.ThreeNode = Ele3Node;
        elseif isempty(strfind(tline,'CPS4')) == 0 % Check elemenent type
            if NDim ~= 2
                error('CPS4 element can only be used in 2D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d '; % Format of 4-Node elements: 5 doubles
            Ele4Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele4Node = Ele4Node{1};
            % Store in structure
            Elements.FourNode = Ele4Node;
        elseif isempty(strfind(tline,'CPS6')) == 0 % Check elemenent type
            if NDim ~= 2
                error('CPS6 element can only be used in 2D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d %d %d ';
            Ele6Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele6Node = Ele6Node{1};
            % Store in structure
            Elements.SixNode = Ele6Node;
        elseif isempty(strfind(tline,'CPS8')) == 0 % Check elemenent type
            if NDim ~= 2
                error('CPS8 element can only be used in 2D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d %d %d %d %d ';
            Ele8Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele8Node = Ele8Node{1};
            
            % Sort nodal connections according to LINKS convention
            Ele8Nodeaux = Ele8Node;
            Ele8Node(:,[3 4 5 6 7 8]) = Ele8Nodeaux(:,[6 3 7 4 8 5]);
            
            % Store in structure
            Elements.EightNode = Ele8Node;
        elseif isempty(strfind(tline,'C3D4')) == 0 % Check elemenent type
            if NDim ~= 3
                error('C3D4 element can only be used in 3D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d ';
            Ele4Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele4Node = Ele4Node{1};
            % Store in structure
            Elements.FourNode = Ele4Node;
        elseif isempty(strfind(tline,'C3D8')) == 0 % Check elemenent type
            if NDim ~= 3
                error('C3D8 element can only be used in 3D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d %d %d %d %d ';
            Ele8Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele8Node = Ele8Node{1};
            % Store in structure
            Elements.EightNode = Ele8Node;
        elseif isempty(strfind(tline,'C3D10')) == 0 % Check elemenent type
            if NDim ~= 3
                error('C3D10 element can only be used in 3D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d %d %d %d %d %d %d ';
            Ele10Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele10Node = Ele10Node{1};
            % Store in structure
            Elements.TenNode = Ele10Node;
        elseif isempty(strfind(tline,'C3D20')) == 0 % Check elemenent type
            if NDim ~= 3
                error('C3D20 element can only be used in 3D analysis')
            end
            % Go to the beginning of first line with matrix elements
            fseek(fileID, filepos, 'bof');

            % Use textscan to read elements
            formatSpec = '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ';
            Ele20Node = textscan(fileID,formatSpec,'Delimiter',',','CollectOutput', true);
            filepos = ftell(fileID);

            % Convert cell array to array of doubles or integers
            Ele20Node = Ele20Node{1};
            
            % The following trick must be done due to the fact that the
            % list of nodes has a line break after column 16
            EleAux = Ele20Node([2:2:length(Ele20Node)] ,:);
            Ele20Node([2:2:length(Ele20Node)] ,:) = [];
            Ele20Node(:,17:21) = EleAux(:,1:5);
            
            % Sort nodal connections according to LINKS convention
            Ele20Nodeaux = Ele20Node;
            Ele20Node(:,[14 15 16 17 18 19 20 21]) = Ele20Nodeaux(:,[18 19 20 21 14 15 16 17]);
            
            % Store in structure
            Elements.TwentyNode = Ele20Node;
        else
            error('Unknown element type found!')
        end
    end
    % Search next type of elements
    tline = fgetl(fileID);
end
end