
%% Load fibre properties database file
FibreMatDat  = 'Datafiles\FibreMaterialData.mat';

fileID = fopen(FibreMatDat,'r');

headerstr = fgets(fileID);
headercell=strsplit(headerstr,' ');

datacell = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f' ,'Delimiter','\n');

fclose(fileID);