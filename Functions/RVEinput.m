%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%            Function that reads the parameters needed to generate        %
%                   a random RVE for hybrid composites.                   %
%                                                                         %
%                 Igor A. Rodrigues Lopes, November 2017                  %
%                              ilopes@fe.up.pt                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       11/2017 - Initial coding                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Input] = RVEinput(RVEDat,fibre1Dat,fibre2Dat)
%% Store data from RVE input file in cell array
disp(' '); disp('Reading RVE parameters')

formatSpec = '%s'; %Read file as string
fileID = fopen(RVEDat,'r');

%Read file and store every line on own line in cell
data = textscan(fileID,formatSpec,'Delimiter','\n');
fclose(fileID);

%% Fibres volume fraction
VolFracStr = strfind(data{1},'Total Fibre Volume Fraction');
for i = 1:length(VolFracStr)
    empty = isempty(cell2mat(VolFracStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.TotalVolFrac = str2double(data{1}(StrLoc));

%% Fibre 2 volume fraction
F2VolStr = strfind(data{1},'Fibre 2 Relative Volume Fraction');
for i = 1:length(F2VolStr)
    empty = isempty(cell2mat(F2VolStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.Fibre2VolFrac = str2double(data{1}(StrLoc));

%% RVE in-plane size
SizeStr = strfind(data{1},'RVE Size');
for i = 1:length(SizeStr)
    empty = isempty(cell2mat(SizeStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.RVEsize = str2double(data{1}(StrLoc));

%% RVE thickness (out-of-plane size)
ThkStr = strfind(data{1},'RVE Thickness');
for i = 1:length(ThkStr)
    empty = isempty(cell2mat(ThkStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.RVEthick = str2double(data{1}(StrLoc));

clear data

%% Get data from fibre 1 input file
fileID = fopen(fibre1Dat,'r');

%Read file and store every line on own line in cell
data = textscan(fileID,formatSpec,'Delimiter','\n');
fclose(fileID);

%% Fibre type 1 radius
RadStr = strfind(data{1},'Fibre Radius');
for i = 1:length(RadStr)
    empty = isempty(cell2mat(RadStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.Radius1 = str2double(data{1}(StrLoc));

clear data

%% Get data from fibre 2 input file
fileID = fopen(fibre2Dat,'r');

%Read file and store every line on own line in cell
data = textscan(fileID,formatSpec,'Delimiter','\n');
fclose(fileID);

%% Fibre type 2 radius
RadStr = strfind(data{1},'Fibre Radius');
for i = 1:length(RadStr)
    empty = isempty(cell2mat(RadStr(i)));
    if empty == 0
        StrLoc = i+1;
    end
end

Input.Radius2 = str2double(data{1}(StrLoc));
end