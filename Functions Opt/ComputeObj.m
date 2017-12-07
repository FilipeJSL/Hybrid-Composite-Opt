function [E] = ComputeObj(MeshFolder, rveFile, s_info)

%% Set path to .out LINKS output file
outFile=strrep(rveFile,'.rve','.out');
outFolder=strrep(rveFile,'.rve','');
OutPath=[MeshFolder '\' outFolder '\' outFile];

%% Data acquisition
delimiterIn = ' ';
headerlinesIn = 1;
data = importdata(OutPath,delimiterIn,headerlinesIn);

%% Save homogenized stress XX (Cauchy) and equivalent strain
SearchKey={'sigma_xx','EquivStrain'}; % Change SearchKey for different targets
Index1=ismember(data.textdata,SearchKey{1});
Index2=ismember(data.textdata,SearchKey{2});

sXX=data.data(:,Index1);
eXX=data.data(:,Index2);

%% Compute the mean elastic modulus
for iincr = 1:s_info.num_incr+1
    E = sXX(iincr)/eXX(iincr);
end

end

