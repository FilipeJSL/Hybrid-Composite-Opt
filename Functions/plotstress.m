function plotstress(Results_Path , irlz)


fileID = fopen(Results_Path,'r');

data = importdata(Results_Path);

fclose(fileID);

%% save homogenised stress_xx (Cauchy) and equivalent strain in vector
sxx = data.data(:,20);
exx = data.data(:,27);

%% plot the stress evolution
figure(1)
plot(exx,sxx)
hold on

%% Compute the mean elastic modulus
ninc = size(sxx);

for iinc = 1:ninc
    figure(2)
    E = sxx(iinc)/exx(iinc);
    plot(irlz,E,'*')
    hold on
end