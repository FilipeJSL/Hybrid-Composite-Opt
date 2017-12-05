function [] = WriteFibreDataFile(Fn,PropStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POTENTIAL UPGRADE: make user option to have "UMATFIBRAND" or "ELASTIC" as
%material types. Here UMATFIBRAND is defined within the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
subfolder='Datafiles\';
strinit='fibre';
strext='.mat';

for i=1:2
    
    fibreFile=[subfolder strinit num2str(i) strext];
    
    fileID=fopen(fibreFile,'w');
    
    FibrePos=PropStruct.Reference==Fn(i);
    
    fprintf(fileID,'** Input file with fibre type %g properties **\n',i);
    fprintf(fileID,'*\n\n');
    fprintf(fileID,'* Using the properties for the fibre number #%g\n\n',PropStruct.Reference(FibrePos));
    fprintf(fileID,'* Fibre Radius - R\n');
    fprintf(fileID,'%f\n\n',PropStruct.Radius(FibrePos));
    fprintf(fileID,'* Material type\n');
    fprintf(fileID,'%s\n\n','UMATFIBRAND'); %Make this a user provided option that comes from the parent script?
    fprintf(fileID,'* Density\n');
    fprintf(fileID,'%g\n\n',0); %Gravitical load option (Not sure if the model is built to take it into account)
    fprintf(fileID,'* List of properties (according to material type)\n');
    fprintf(fileID,'74E+3  0.2 (not used for UMAT)\n\n');
    fprintf(fileID,'* Number of Hardening points\n');
    fprintf(fileID,'%g\n\n',12); %Do not change this value. It's associated with the no. of properties defined in UMATFRIBRAND material type
    fprintf(fileID,'* List of Hardening points\n');
    fprintf(fileID,'%-11g# Longitudinal Young modulus E11 [MPa]\n',PropStruct.LYoungMod(FibrePos));
    fprintf(fileID,'%-11g# Transverse Young modulus E22 [MPa]\n',PropStruct.TYoungMod(FibrePos));
    fprintf(fileID,'%-11g# Major Poisson coefficient nu12\n',PropStruct.PoissonCoef(FibrePos));
    fprintf(fileID,'%-11g# Longitudinal shear modulus G12 [MPa]\n',PropStruct.LShearMod(FibrePos));
    fprintf(fileID,'%-11g# Transverse shear modulus G23 [MPa]\n',PropStruct.TShearMod(FibrePos));
    fprintf(fileID,'%-11g# Initial tension failure stress sigma0 [MPa]\n',PropStruct.Sigma0(FibrePos));
    fprintf(fileID,'%-11g# Longitudinal thermal expansion alpha11\n',PropStruct.LThermalExp(FibrePos));
    fprintf(fileID,'%-11g# Transverse thermal expansion alpha22\n',PropStruct.TThermalExp(FibrePos));
    fprintf(fileID,'%-11g# Fracture toughness Gcf [N/mm]\n',PropStruct.FracThoughness(FibrePos));
    fprintf(fileID,'%-11g# Weibull parameter m\n',PropStruct.WeibullPar(FibrePos));
    fprintf(fileID,'%-11g# Characteristic fibre length [mm]\n',PropStruct.CharacFibreLeng(FibrePos));
    fprintf(fileID,'           # Actual fibre length [mm]\n'); %This property is written later on in a subsequent function (it depends on the RVE thickness)

    fclose(fileID);
end

end

