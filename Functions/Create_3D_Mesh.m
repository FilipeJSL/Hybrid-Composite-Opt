%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Function that writes Python script                     %
%                   to generate a 3D mesh in Abaqus                       %
%                                                                         %
%                Created by Igor Lopes October 2017 ©                     %
%                                                                         %
%  Based on functions "f_mesh_quad_per_3D.m" and "Create_2D_Mesh"         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Changelog:                                                            %
%                                                                         %
%       10/2017 - First Version                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [status_mesh] = Create_3D_Mesh(dir_name)

global a b c Fibre_pos N_fibre S_base Rmax Abaqus;

%% Create Python File for Abaqus CAE with Fibre Distribution

disp('Started creation of Python file');

file_name = strcat(dir_name,'/mesh_3D.py');
fid = fopen(file_name,'wt');

%% Heading

fprintf(fid,'# -*- coding: mbcs -*-\n');
fprintf(fid,'#\n');
fprintf(fid,'# Abaqus/CAE \n');
fprintf(fid,'# Created on %s\n',datestr(now));
fprintf(fid,'#\n');
fprintf(fid,'\n');
fprintf(fid,'# from driverUtils import executeOnCaeGraphicsStartup\n');
fprintf(fid,'# executeOnCaeGraphicsStartup()\n');
fprintf(fid,'#: Executing "onCaeGraphicsStartup()" in the site directory ...\n');
fprintf(fid,'from abaqus import *\n');
fprintf(fid,'from abaqusConstants import *\n');
fprintf(fid,'session.Viewport(name=''Viewport: 1'', origin=(0.0, 0.0), width=196.9, \n');
fprintf(fid,'    height=171.7)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].makeCurrent()\n');
fprintf(fid,'session.viewports[''Viewport: 1''].maximize()\n');
fprintf(fid,'from caeModules import *\n');
fprintf(fid,'from driverUtils import executeOnCaeStartup\n');
fprintf(fid,'executeOnCaeStartup()\n');
fprintf(fid,'Mdb()\n');

%% Fibre Sketches

fprintf(fid,'#\n');
fprintf(fid,'# Fibres Sketches\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%18.15f)\n',a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
for i=1:1:N_fibre
    XC = Fibre_pos(i,2); YC = Fibre_pos(i,3);
    fprintf(fid,'s.CircleByCenterPerimeter(center=(%18.15f,%18.15f), point1=(\n',XC,YC);
    fprintf(fid,'    %18.15f, %18.15f))\n',XC+Fibre_pos(i,6),YC);
end
fprintf(fid,'mdb.models[''Model-1''].ConstrainedSketch(name=''Fibre_Sketch_Prov'', objectToCopy=s)\n');
fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
fprintf(fid,'    toName=''Fibre_Sketch_Prov'')\n');
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.rectangle(point1=(0.0, 0.0), point2=(%7.5f, %7.5f))\n',a,b);
fprintf(fid,'s.rectangle(point1=(%7.5f, %7.5f), point2=(%7.5f, %7.5f))\n',-2*Rmax,-2*Rmax,a+2*Rmax,b+2*Rmax);
fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
fprintf(fid,'    toName=''Fibre_Sketch_Trim'')\n');
fprintf(fid,'s.unsetPrimaryObject()\n');

%% Matrix Sketch

fprintf(fid,'#\n');
fprintf(fid,'# Matrix Sketch\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.rectangle(point1=(0.0, 0.0), point2=(%7.5f, %7.5f))\n',a,b);
fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
fprintf(fid,'    toName=''Matrix_Sketch'')\n');
fprintf(fid,'s.unsetPrimaryObject()\n');

%% Provisional Fibre Part

fprintf(fid,'#\n');
fprintf(fid,'# Fibres Part (Provisional)\n');
fprintf(fid,'#\n');
fprintf(fid,'s1 = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',a);
fprintf(fid,'g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints\n');
fprintf(fid,'s1.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s1.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
fprintf(fid,'s1.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Fibre_Sketch_Prov''])\n');
fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Fibre_Part_Prov'', dimensionality=THREE_D, \n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'p.BaseSolidExtrude(sketch=s1, depth=%7.5f)\n',c);
fprintf(fid,'s1.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''Fibre_Sketch_Prov'']\n');

%% Fibre Part Trimmer Only

fprintf(fid,'#\n');
fprintf(fid,'# Fibres Trimmer Part\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
fprintf(fid,'s.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Fibre_Sketch_Trim''])\n');
fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Fibre_Part_Trim'', dimensionality=THREE_D, \n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'p.BaseSolidExtrude(sketch=s, depth=%7.5f)\n',c);
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''Fibre_Sketch_Trim'']\n');

%% Matrix Part in Cube Shape Only

fprintf(fid,'#\n');
fprintf(fid,'# Matrix Part (Cube only)\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
fprintf(fid,'s.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Matrix_Sketch''])\n');
fprintf(fid,'session.viewports[''Viewport: 1''].view.fitView()\n');
fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Matrix_Part'', dimensionality=THREE_D, \n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
fprintf(fid,'p.BaseSolidExtrude(sketch=s, depth=%7.5f)\n',c);
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=p)\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''Matrix_Sketch'']\n');

%% Creation of Real Fibre Part

fprintf(fid,'#\n');
fprintf(fid,'# Fibre Part\n');
fprintf(fid,'#\n');
fprintf(fid,'a = mdb.models[''Model-1''].rootAssembly\n');
fprintf(fid,'a.DatumCsysByDefault(CARTESIAN)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'a.Instance(name=''Fibre_Instance_Prov'', part=p, dependent=OFF)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'a.Instance(name=''Fibre_Instance_Trim'', part=p, dependent=OFF)\n');
fprintf(fid,'a.InstanceFromBooleanCut(name=''Fibre_Final_Part'', \n');
fprintf(fid,'    instanceToBeCut=mdb.models[''Model-1''].rootAssembly.instances[''Fibre_Instance_Prov''], \n');
fprintf(fid,'    cuttingInstances=(a.instances[''Fibre_Instance_Trim''], ), originalInstances=DELETE)\n');
fprintf(fid,'a.makeIndependent(instances=(a.instances[''Fibre_Final_Part-1''], ))\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.features.changeKey(\n');
fprintf(fid,'    fromName=''Fibre_Final_Part-1'', toName=''Fibre_Final_Instance'')\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');

%% Creation of Real Matrix Part

fprintf(fid,'#\n');
fprintf(fid,'# Matrix Part\n');
fprintf(fid,'#\n');
fprintf(fid,'a = mdb.models[''Model-1''].rootAssembly\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
fprintf(fid,'a.Instance(name=''Matrix_Instance'', part=p, dependent=OFF)\n');
fprintf(fid,'a.InstanceFromBooleanCut(name=''Matrix_Final_Part'', \n');
fprintf(fid,'    instanceToBeCut=mdb.models[''Model-1''].rootAssembly.instances[''Matrix_Instance''], \n');
fprintf(fid,'    cuttingInstances=(a.instances[''Fibre_Final_Instance''], ), originalInstances=DELETE)\n');
fprintf(fid,'a.makeIndependent(instances=(a.instances[''Matrix_Final_Part-1''], ))\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.features.changeKey(\n');
fprintf(fid,'    fromName=''Matrix_Final_Part-1'', toName=''Matrix_Final_Instance'')\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Matrix_Part'']\n');

%% Creation of Final Geometry

fprintf(fid,'#\n');
fprintf(fid,'# Assembly and Meshing\n');
fprintf(fid,'#\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Final_Part'']\n');
fprintf(fid,'a.Instance(name=''Fibre_Final_Instance'', part=p, dependent=OFF)\n');
fprintf(fid,'a.InstanceFromBooleanMerge(name=''Final_Stuff'', instances=(\n');
fprintf(fid,'    a.instances[''Matrix_Final_Instance''], a.instances[''Fibre_Final_Instance''], ), \n');
fprintf(fid,'    keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Matrix_Final_Part'']\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Final_Part'']\n');
fprintf(fid,'a.makeIndependent(instances=(a.instances[''Final_Stuff-1''], ))\n');
fprintf(fid,'partInstances =(a.instances[''Final_Stuff-1''], )\n');
fprintf(fid,'a.seedPartInstance(regions=partInstances, size=%7.5f, deviationFactor=0.1)\n',S_base);
fprintf(fid,'elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)\n');
fprintf(fid,'elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)\n');
fprintf(fid,'elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)\n');
fprintf(fid,'c1 = a.instances[''Final_Stuff-1''].cells\n');
fprintf(fid,'pickedRegions =(c1, )\n');
fprintf(fid,'a.setElementType(regions=pickedRegions, elemTypes=(elemType1,))\n');
fprintf(fid,'allCells = a.instances[''Final_Stuff-1''].cells\n');
fprintf(fid,'a.setMeshControls(regions=allCells, technique=SWEEP, allowMapped=False, \n');
fprintf(fid,'    elemShape=HEX)\n');
fprintf(fid,'a.generateMesh(regions=partInstances, seedConstraintOverride=OFF,\n');
fprintf(fid,'    meshTechniqueOverride=OFF)\n');

%% View settings

fprintf(fid,'#\n');
fprintf(fid,'# Viewports\n');
fprintf(fid,'#\n');
fprintf(fid,'session.viewports[''Viewport: 1''].assemblyDisplay.setValues(mesh=ON)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].view.fitView()\n');
fprintf(fid,'session.viewports[''Viewport: 1''].view.setProjection(projection=PARALLEL)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=a)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].assemblyDisplay.meshOptions.setValues(\n');
fprintf(fid,'    meshTechnique=ON)\n');

%% Export mesh to .inp file

fprintf(fid,'#\n');
fprintf(fid,'# Job creation\n');
fprintf(fid,'#\n');
fprintf(fid,'mdb.Job(name=''mesh'', model=''Model-1'', type=ANALYSIS, explicitPrecision=SINGLE, \n');
fprintf(fid,'    nodalOutputPrecision=SINGLE, description='''', \n');
fprintf(fid,'    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, \n');
fprintf(fid,'    numDomains=1, userSubroutine='''', numCpus=1, memory=90, \n');
fprintf(fid,'    memoryUnits=PERCENTAGE, scratch='''', echoPrint=OFF, modelPrint=OFF, \n');
fprintf(fid,'    contactPrint=OFF, historyPrint=OFF)\n');
fprintf(fid,'import os\n');
fprintf(fid,'os.chdir(r''%s'')\n',dir_name);
fprintf(fid,'mdb.jobs[''mesh''].writeInput(consistencyChecking=OFF)\n');

fclose(fid);

disp(' ');
disp('Creation of Python file COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);

%% Execute Python Script in ABAQUS CAE

disp('Started generation of mesh'); disp(' ');

dos(strcat(Abaqus,' cae noGUI=',dir_name,'\mesh_3D.py'));
disp(' ');
disp('Generation of mesh COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);
status_mesh = 1;

% Delete Abaqus output files
delete('*.rpy')
delete('*.log')
delete('*.rec')
end