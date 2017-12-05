function [PropStruct] = GetFibreProps( FibreMatDat )

fileID = fopen(FibreMatDat,'r');

headerstr = fgetl(fileID);
headercell=strsplit(headerstr,' ');

formatSpec='%f %f %f %f %f %f %f %f %f %f %f %f %f %f';

datacell = textscan(fileID,formatSpec,'Delimiter','\n');

fclose(fileID);

%Stores the fibre properties database as a structure
PropStruct=cell2struct(datacell,headercell,2);

end

