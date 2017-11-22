%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Check for Incompatibility of the Fibre Position       %
%                                                                         %
%  This file is part of the RAND_uSTRU_GEN.m file                         %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v2.0                                   %
%                                                                         %
%                Ant�nio Rui Melro - antonio.melro@fe.up.pt               %
%                                June 2006                                %
%                                                                         %
%		Rodrigo Paiva Tavares - rptavares@inegi.up.pt		  %
%				October 2015				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [fibre_overlap,min] = f_overlap(X_TMP,Y_TMP,fibre2test,r,ftype)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global N_fibre Fibre_pos DISTMIN a b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify Distribution for Fibre Overlaps                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fibre_overlap = 0;
min = a*2;
for k=1:N_fibre
    DISTMIN_2=Fibre_pos(k,6)+r+(DISTMIN-2)*(Fibre_pos(k,6)+r)/2;
xlim0 = X_TMP - 4*DISTMIN_2; xlim1 = X_TMP + 4*DISTMIN_2;
ylim0 = Y_TMP - 4*DISTMIN_2; ylim1 = Y_TMP + 4*DISTMIN_2;
    if k ~= fibre2test
        xx = Fibre_pos(k,2); yy = Fibre_pos(k,3);
        if xx > xlim0 && xx < xlim1 && yy > ylim0 && yy < ylim1
            new_dist = sqrt((xx-X_TMP)^2 + (yy-Y_TMP)^2);
            if new_dist < DISTMIN_2
                fibre_overlap = 1;
            else
                if new_dist < min
                    min = new_dist;
                end
            end
        end
    end
    
end

%  if ftype==0
%      if Y_TMP >3/4*b+r || Y_TMP< 1/4*b -r 
%          fibre_overlap=1;
%      end
%  else
%      if Y_TMP <3/4*b-r && Y_TMP> 1/4*b+r  
%          fibre_overlap=1;
%      end
%  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

