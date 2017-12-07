function [rveFile] = RunLINKS(MeshFolder, FolderName, CurrentFolder, s_info)

%% Write batch file to run LINKS
batName='run_LINKS.bat';
BatPath=[FolderName '\' batName];
fileID=fopen(BatPath,'w');

Exe=['"' CurrentFolder '\LINKS.exe"'];
rveFile=['RVE_' s_info.BC{1} '.rve'];
InputPath=['"' CurrentFolder '\' MeshFolder  '\' rveFile '"']; %Is it possible that several BC are stated, instead of just one? (See LINKS_DatFile.m line 56)

formatSpec='%s %s';

fprintf(fileID,formatSpec,Exe,InputPath);

fclose(fileID);

%% Call to run_LINKS.bat
[~,~]=dos([CurrentFolder '\' FolderName '\' batName]); %The [~,~] prevents Matlab from printing in the command window

end

