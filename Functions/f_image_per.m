%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Generate an image of the RVE in .eps format           %
%                                                                         %
%  This file is part of the XXXX_PER_uSTRU_GEN.m sequence of files        %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v2.0                                   %
%                                                                         %
%                António Rui Melro - antonio.melro@fe.up.pt              %
%                             October  2008                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   MODIFICATION LOG                                                      %
%                                                                         %
%         10/2016 - Modified by Lars Peperkamp.                           %
%                   No longer show image, but save final image as .pdf.   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Legend for variable "index1":                                          %
%       0 - image in the end of routine RAND_PER_uSTRU_GEN.m              %
%       1 - image in the end of step 1 of routine RAND_PER_uSTRU_GEN.m    %
%       2 - image in the end of step 2 of routine RAND_PER_uSTRU_GEN.m    %
%       3 - image in the end of step 3 of routine RAND_PER_uSTRU_GEN.m    %
%       4 - image for SPER_uSTRU_GEN.m                                    %
%       5 - image for HPER_uSTRU_GEN.m                                    %
%       6 - image for WONG_uSTRU_GEN.m                                    %
%       9 - image for animation sequence of routine RAND_PER_uSTRU_GEN.m  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function f_image_per(index1,FolderName)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%dir_name%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global a b N_cycles Fibre_pos N_fibre vorim_option Vol_fibre A_1_fibre...
    A_2_fibre Rmax Rmin;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Opens Figure and Adjusts Properties                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure('visible','on');
hold on;
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'TickDir','out');
set(gca,'Units','points');
rectangle('Position',[0,0,a,b],'EdgeColor','w');
axis([-2*Rmax a+2*Rmax -2*Rmax b+2*Rmax]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the Title                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if index1 == 1 || index1 == 2 || index1 == 3
    title1 = num2str(N_cycles);
    title2 = num2str(index1);
    title3 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Step: ', title2,';   v_f=', title3,'%']);
elseif index1 == 0
    title1 = num2str(N_cycles-1);
    title3 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Step: End;   v_f=', title3,'%']);
elseif index1 == 4 || index1 == 5
    title3 = num2str(Vol_fibre*100);
    title(['Fibre Volume = ', title3,'%']);
elseif index1 == 6
    title1 = num2str(N_cycles);
    title2 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Fibre Volume = ', title2,'%']);
elseif index1 == 9
    title3 = num2str(Vol_fibre*100);
    title(['v_f=', title3,'%']);
end
whitebg('black');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the Position of all Fibres                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Fibre_scale = 72^2*25.4^2*(Rmin/Rmax);

for i=1:N_fibre
    if Fibre_pos(i,5)==0
        scatter(Fibre_pos(i,2),Fibre_pos(i,3),A_1_fibre*Fibre_scale,...
            [0.8,0.8,0.8],'filled');
    end
    if Fibre_pos(i,5)==1
        scatter(Fibre_pos(i,2),Fibre_pos(i,3),A_2_fibre*Fibre_scale,...
            [0.4,0.4,0.4],'filled');
    end
end
if vorim_option == 1
    voronoi(Fibre_pos(:,2),Fibre_pos(:,3),'w-');
end
hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves the Image in .eps Format                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if index1 == 0
    %DirName = strcat(cd,'\Final_Fibre_Distribution_Image\');
    %if exist(DirName,'dir') == 0, mkdir(DirName); end
    %name = strcat(DirName,FolderName);
    %name = name{1};
    name = strcat(FolderName,'\Distribution');
    saveas(1,name,'pdf');
    %     print -dps2 -append anim;
    close(1);
end
if index1 == 1
    name = num2str(3*N_cycles-2,'%03d');
    saveas(1,name,'epsc');
    %     print -dps2 -append anim;
    close(1);
end
if index1 == 2
    name = num2str(3*N_cycles-1,'%03d');
    saveas(1,name,'epsc');
    %     print -dps2 -append anim;
    close(1);
end
if index1 == 3
    name = num2str(3*N_cycles,'%03d');
    saveas(1,name,'epsc');
    %     print -dps2 -append anim;
    close(1);
end
if index1 == 4
    name = 'SPER_uSTRU';
    saveas(1,name,'epsc');
    close(1);
end
if index1 == 5
    name = 'HPER_uSTRU';
    saveas(1,name,'epsc');
    close(1);
end
if index1 == 6
    name = 'WONG_uSTRU';
    saveas(1,name,'epsc');
    close(1);
end
if index1 == 9
    print -dps2 -append animation;
    close(1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%