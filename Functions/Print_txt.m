%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              Function that prints information in .txt file              %
%                                                                         %
%                Created by Lars Peperkamp October 2016 ©                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       10/2016 - First Version                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Print_txt(FilePath,FileName,Input,Format)
%% Combine FileName and FilePath to FullFilePath
FullFilePath = strcat(FilePath,FileName);

%% Write Input in defined Format
fileID = fopen(FullFilePath,'w');
for i = 1:length(Input)
    fprintf(fileID,Format,Input(i,:));
end
fclose(fileID);

end