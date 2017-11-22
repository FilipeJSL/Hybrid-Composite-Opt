%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Generate a Random Distribution of Circles             %
%                                                                         %
%  This file is part of the FEA_Automatic.m sequence of files             %
%  and should not be used separately from it.                             %
%                                                                         %
%                                 v12.0                                   %
%                                                                         %
%                Antonio Rui Melro - antonio.melro@fe.up.pt               %
%                             October  2008                               %  
%               Rodrigo Paiva Tavares - em10140@fe.up.pt
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   MODIFICATION LOG                                                      %
%         
%         2014/11/** - Modification to allow hybrid composites   
%         2008/10/25 - Reorganised all process of calling scripts         %
%         2007/12/09 - Removed control/check variables                    %
%         2007/03/28 - Ammended STEP four                                 %
%         2006/11/17 - Added possibility to generate mesh                 %
%         2006/07/13 - Divided the routine in different functions         %
%         2006/07/04 - Optimized the output process by including arcs     %
%         2006/07/01 - Added output instructions to generate Python file  %
%         2006/06/30 - Added calculus of std_dist_adim of Voronoi Polygon %
%                      Distances                                          %
%         2006/06/29 - Added calculus of std_area_adim of Voronoi Polygon %
%                      Areas and parametrized all routine                 %
%         2006/06/26 - Added possibility to analyse Interlaminar Area     %
%         2006/06/20 - Removed First Heuristic                            %
%         2006/06/20 - Optimisation of Third Heuristic                    %
%         2006/06/19 - Optimisation of All Heuristics                     %
%         2006/06/17 - Implementation of Third Heuristic                  %
%         2006/06/16 - Optimisation of First and Second Heuristics        %
%         2006/06/15 - Implementation of Second Heuristic                 %
%         2006/06/15 - Optimisation of First Heuristic                    %
%         2006/06/14 - Implementation of First Heuristic                  %
%         2006/06/13 - First version
%        
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   APPENDED FUNCTIONS                                                    %
%         f_overlap - Function to check for the existence of an overlap   %
%                  between the fibre of interest and the fibres in an     %
%                  area around the fibre of interest                      %
%         f_image_per - Generates image files of fibre distribution       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [status]=RAND_PER_uSTRU_GEN_3D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global R Vol_fibre_req DISTMIN N_cycles Max_fibres N_guesses_max ...
    N_cycles_max N_change Square_size Square_inc inter_option ...
    image_option a b A_1_fibre Vol_fibre Fibre_pos N_fibre S_base...
    Fibre_type_1 R1 A_2_fibre Rmax Rmin N_fibre_real;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of Auxiliary Parameters                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if inter_option == 1
    A_total = a*(b-2.0*R); % total area of image minus the interlaminar area
else
    A_total = a*b;         % total area of image
end
A_1_fibre = pi*R^2;        % area of a single fibre type 1
A_2_fibre = pi*R1^2;       % area of a single fibre type 2
a0 = -1*Rmin; b0 = a0;        % initial coordinate for position of fibre
af = a + Rmin; bf = b + Rmin;    % final coordinate for position of fibre
N_cycles = 1;              % number of cycles of the algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Output Variable                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fibre_pos  % Matrix of coordinates
%                      % 1st Column = fibre ID
%                      % 2nd Column = X-coordinate of centre
%                      % 3rd Column = Y-coordinate of centre
%                      % 4th Column = Splitted fibre:
%                      0=No;2=half;4=corners
%                      % 5th Column = Fibre Type     
%                      % 6th Column =  Fibre radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of Calculus                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Vol_fibre = 0;  Vol_fibre_1=0;    Vol_fibre_2=0;     N_fibre_real = 1;
Vec_mem = zeros(Max_fibres,N_cycles_max+1);  % Memory vector for First Heuristic
N_fibre = 0;N_fibre_1=0;N_fibre_2=0;         % number of fibres
N_fibre = N_fibre + 1;                       % Updated number of fibres
%
if Fibre_type_1==1
    ft=0;
elseif Fibre_type_1==0
    ft=1;
else
    if R>R1
        ft=0;
    else
        ft=1;
    end
end
Fibre_pos(1,5)=ft;
if Fibre_pos(1,5)==0
    Fibre_pos(1,6)=R;
    Vol_fibre = Vol_fibre + A_1_fibre / A_total;
    Vol_fibre_1=Vol_fibre_1+A_1_fibre / A_total;
elseif Fibre_pos(1,5)==1
    Fibre_pos(1,6)=R1;
    Vol_fibre = Vol_fibre + A_2_fibre / A_total;
    Vol_fibre_2=Vol_fibre_2+A_2_fibre / A_total;
end
XC = unifrnd(Fibre_pos(1,6),a-Fibre_pos(1,6));  % X-position of first fibre
YC = unifrnd(Fibre_pos(1,6),b-Fibre_pos(1,6));  % Y-position of first fibre
Fibre_pos(1,2) = XC; Fibre_pos(1,3) = YC; Fibre_pos(1,4) = 0;

N_attempts = 0;                % number of attempts on fibre collocation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomly Determine Position of Fibres and Compatibility Check           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
while Vol_fibre_req > Vol_fibre
    disp(' '); disp('Iteration number: '); disp(N_cycles);
    while N_attempts < N_cycles*N_guesses_max && Vol_fibre_req > Vol_fibre && N_fibre < Max_fibres
        N_attempts = N_attempts + 1; split = 0;
        %XC = unifrnd(a0,af); YC = unifrnd(b0,bf);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Fibre_type_1==1
            ft=0;
        elseif Fibre_type_1==0
            ft=1;
        elseif Vol_fibre_1< Fibre_type_1*(Vol_fibre)
            ft=0;     
        elseif Vol_fibre_2< (1-Fibre_type_1)*(Vol_fibre)
            ft=1;
        end
        if ft==0
            r=R;
            A_fibre=A_1_fibre;                  
        else
            r=R1;
            A_fibre=A_2_fibre;
        end
        
        
        XC = unifrnd(-r,a+r); YC = unifrnd(-r,b+r);

        
%         if ft==0
%             XC = unifrnd(-r,a+r); YC = unifrnd(b/4-r,3/4*b+r);
%         else
%             XC = unifrnd(-r,a+r);
%             if rand>0.5
%                 YC = unifrnd(-r,1/4*b+r);
%             else
%                 YC = unifrnd(3/4*b-r,b +r);
%             end
%          end
                %%%%%%%%%%%%%%%%%%%%%%%%%
        if XC >= r && XC <= a-r && YC >= r && YC <= b-r
            [check dist] = f_overlap(XC,YC,N_fibre + 1,r,ft);
            if check == 0
                N_fibre=N_fibre+1;
                %%%%%%%%%%%%%%%%%%%%%%%
                Fibre_pos(N_fibre,5)=ft;
                Fibre_pos(N_fibre,6)=r;
                if Fibre_pos(N_fibre,5)==0
                    N_fibre_1=N_fibre_1+1;
                    Vol_fibre_1=Vol_fibre_1+A_fibre/A_total;
                elseif Fibre_pos(N_fibre,5)==1
                    N_fibre_2=N_fibre_2+1;
                    Vol_fibre_2=Vol_fibre_2+A_fibre/A_total;
                end
                %%%%%%%%%%%%%%%%%%%
                Fibre_pos(N_fibre,2) = XC; Fibre_pos(N_fibre,3) = YC;
                Fibre_pos(N_fibre,4) = 0;
              
                Vol_fibre = Vol_fibre + A_fibre / A_total;
                N_fibre_real = N_fibre_real + 1;
            end
        elseif XC >= r && XC <= a-r && YC < r    %Z2
            XC1 = XC; YC1 = YC + b; split = 2;
        elseif XC >= r && XC <= a-r && YC > b-r  %Z3
            XC1 = XC; YC1 = YC - b; split = 2;
        elseif YC >= r && YC <= b-r && XC < r    %Z4
            XC1 = XC + a; YC1 = YC; split = 2;
        elseif YC >= r && YC <= b-r && XC > a-r  %Z5
            XC1 = XC - a; YC1 = YC; split = 2;
        elseif XC < r && YC < r                  %Z6
            XC1 = XC + a; YC1 = YC + b; split = 4;
        elseif XC > a-r && YC < r                %Z7
            XC1 = XC - a; YC1 = YC + b; split = 4;
        elseif XC > a-r && YC > b-r              %Z8
            XC1 = XC - a; YC1 = YC - b; split = 4;
        elseif XC < r && YC > b-r                %Z9
            XC1 = XC + a; YC1 = YC - b; split = 4;
        end
        %
        if split == 2
            [check dist] = f_overlap(XC,YC,N_fibre+1,r,ft);
            if check == 0
                [check dist] = f_overlap(XC1,YC1,N_fibre+1,r,ft);
                if check == 0
                    N_fibre = N_fibre + 1;
                    Fibre_pos(N_fibre,2) = XC; Fibre_pos(N_fibre,3) = YC;
                    Fibre_pos(N_fibre,4) = 2;
                    %%%%%%%%%%%%%%%%%%%%%%%
                    Fibre_pos(N_fibre,5)=ft;
                    Fibre_pos(N_fibre,6)=r;
                    if Fibre_pos(N_fibre,5)==0
                        N_fibre_1=N_fibre_1+1;
                        Vol_fibre_1=Vol_fibre_1+A_fibre/A_total;
                    elseif Fibre_pos(N_fibre,5)==1
                        N_fibre_2=N_fibre_2+1;
                        Vol_fibre_2=Vol_fibre_2+A_fibre/A_total;
                    end
                %%%%%%%%%%%%%%%%%%%
                    N_fibre = N_fibre + 1;
                    Fibre_pos(N_fibre,2) = XC1; Fibre_pos(N_fibre,3) = YC1;
                    Fibre_pos(N_fibre,4) = 2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    Fibre_pos(N_fibre,5)=ft;
                    Fibre_pos(N_fibre,6)=r;
                    %%%%%%%%%%%%%%%%
                    Vol_fibre = Vol_fibre + A_fibre / A_total;
                    N_fibre_real = N_fibre_real + 1;
                end
            end
        elseif split == 4
            [check dist] = f_overlap(XC,YC,N_fibre + 1,r,ft);
            if check == 0
                [check dist] = f_overlap(XC,YC1,N_fibre + 1,r,ft);
                if check == 0
                    [check dist] = f_overlap(XC1,YC,N_fibre + 1,r,ft);
                    if check == 0
                        [check dist] = f_overlap(XC1,YC1,N_fibre + 1,r,ft);
                        if check == 0
                            N_fibre = N_fibre + 1;
                            Fibre_pos(N_fibre,2) = XC; Fibre_pos(N_fibre,3) = YC;
                            Fibre_pos(N_fibre,4) = 4;
                            %%%%%%%%%%%%%%%%%%%%%%%
                            Fibre_pos(N_fibre,5)=ft;
                            Fibre_pos(N_fibre,6)=r;
                            if Fibre_pos(N_fibre,5)==0
                                N_fibre_1=N_fibre_1+1;
                                Vol_fibre1_=Vol_fibre_1+A_fibre/A_total;
                            elseif Fibre_pos(N_fibre,5)==1
                                N_fibre_2=N_fibre_2+1;
                                Vol_fibre_2=Vol_fibre_2+A_fibre/A_total;
                            end
                            %%%%%%%%%%%%%%%%%%%
                            N_fibre = N_fibre + 1;
                            Fibre_pos(N_fibre,2) = XC; Fibre_pos(N_fibre,3) = YC1;
                            Fibre_pos(N_fibre,4) = 4;
                            %%%%%%%%%%%%%%%%%%%%%%%%%
                            Fibre_pos(N_fibre,5)=ft;
                            Fibre_pos(N_fibre,6)=r;
                            %%%%%%%%%%%%%%%%
                            N_fibre = N_fibre + 1;
                            Fibre_pos(N_fibre,2) = XC1; Fibre_pos(N_fibre,3) = YC;
                            Fibre_pos(N_fibre,4) = 4;
                            %%%%%%%%%%%%%%%%%%%%%%%%%
                            Fibre_pos(N_fibre,5)=ft;
                            Fibre_pos(N_fibre,6)=r;
                            %%%%%%%%%%%%%%%%
                            N_fibre = N_fibre + 1;
                            Fibre_pos(N_fibre,2) = XC1; Fibre_pos(N_fibre,3) = YC1;
                            Fibre_pos(N_fibre,4) = 4;
                               %%%%%%%%%%%%%%%%%%%%%%%%%
                            Fibre_pos(N_fibre,5)=ft;
                            Fibre_pos(N_fibre,6)=r;
                                %%%%%%%%%%%%%%%%
                            Vol_fibre = Vol_fibre + A_fibre / A_total;
                            N_fibre_real = N_fibre_real + 1;
                        end
                    end
                end
            end
        end
    end
    disp('Time elapsed after step 1 [min]: '); disp(toc/60);
    format('long');
    disp('Fibre Volume achieved after step 1: '); disp(Vol_fibre);
    disp('Achieved Type 1 Fibre Volume: '); disp(Vol_fibre_1);
    disp('Achieved Type 2 Fibre Volume: '); disp(Vol_fibre_2);
    format('short');
        if image_option == 1
            index1 = 1;
            f_image_per(index1);
        end
    %
    % First Heuristic - Move Fibres Around to Gain More Empty Areas
    %
    if N_fibre <= Max_fibres
        fibre_jump = 0;
        for i=1:N_fibre
            if fibre_jump ~= 0
                fibre_jump = fibre_jump - 1; continue
            end
            if i>N_fibre, break, end
            X_TMP = Fibre_pos(i,2); Y_TMP = Fibre_pos(i,3);
            MIN = a*2; IMIN = 0;
            for j=1:N_fibre
                if i ~= j
                    x_delta = abs(Fibre_pos(j,2)-X_TMP);
                    y_delta = abs(Fibre_pos(j,3)-Y_TMP);
                    if x_delta < 4*(Fibre_pos(i,6)+Fibre_pos(j,6)+(DISTMIN-2)*(Fibre_pos(i,6)+Fibre_pos(j,6))/2) && y_delta < 4*(Fibre_pos(i,6)+Fibre_pos(j,6)+(DISTMIN-2)*(Fibre_pos(i,6)+Fibre_pos(j,6))/2)
                        DIST_TMP = sqrt(x_delta^2 + y_delta^2);
                        if DIST_TMP < MIN
                            if N_cycles == 1
                                IMIN = j;
                                MIN = DIST_TMP;
                                Vec_mem(i,N_cycles+1) = j;
                            elseif mod(N_cycles,N_change) == 0
                                if Vec_mem(i,N_cycles) ~= j || Vec_mem(i,N_cycles-1) ~= j
                                    IMIN = j;
                                    MIN = DIST_TMP;
                                    Vec_mem(i,N_cycles+1) = j;
                                end
                            else
                                if Vec_mem(i,N_cycles) ~= j
                                    IMIN = j;
                                    MIN = DIST_TMP;
                                    Vec_mem(i,N_cycles+1) = j;
                                end
                            end
                        end
                    end
                end
            end
            Delta = unifrnd(0,1);
            if IMIN ~= 0
                x_delta = Fibre_pos(IMIN,2)-X_TMP; y_delta = Fibre_pos(IMIN,3)-Y_TMP;
                DISTMIN_2=Fibre_pos(i,6)+Fibre_pos(IMIN,6)+(DISTMIN-2)*(Fibre_pos(i,6)+Fibre_pos(IMIN,6))/2;
            end
            delta_len = sqrt(x_delta^2 + y_delta^2); k = 1 - DISTMIN_2/delta_len;
            X_TMP = Fibre_pos(i,2) + Delta*k*x_delta;
            Y_TMP = Fibre_pos(i,3) + Delta*k*y_delta;
            [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
            if check == 0
                if Fibre_pos(i,4) == 0     % Z1
                    if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;   %Z1->Z1
                    else
                        split = 0;
                        if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)        % Z1->Z2
                            X_TMP1 = X_TMP; Y_TMP1 = Y_TMP + b; split = 2;
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)  % Z1->Z3
                            X_TMP1 = X_TMP; Y_TMP1 = Y_TMP - b; split = 2;
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)    % Z1->Z4
                            X_TMP1 = X_TMP + a; Y_TMP1 = Y_TMP; split = 2;
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP > a-Fibre_pos(i,6)  % Z1->Z5
                            X_TMP1 = X_TMP - a; Y_TMP1 = Y_TMP; split = 2;
                        end
                        if split == 2
                            [check dist] = f_overlap(X_TMP1,Y_TMP1,i,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 2;
                                N_fibre = N_fibre + 1;
                                for j=N_fibre:-1:i+2
                                    for jj=2:6, Fibre_pos(j,jj) = Fibre_pos(j-1,jj); end
                                    for jj=1:N_cycles_max+1, Vec_mem(j,jj) = Vec_mem(j-1,jj); end
                                end
                                Fibre_pos(i+1,2) = X_TMP1; Fibre_pos(i+1,3) = Y_TMP1;
                                Fibre_pos(i+1,4) = 2; Fibre_pos(i+1,5)=Fibre_pos(i,5); Fibre_pos(i+1,6)=Fibre_pos(i,6);
                            end
                        end
                    end
                elseif Fibre_pos(i,4) == 2
                    if i ~= N_fibre
                        if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                            ics = 1;
                        else
                            ics = -1;
                        end
                    else ics = -1;
                    end
                    if XC >= Fibre_pos(i,6) && XC <= a-Fibre_pos(i,6) && YC < Fibre_pos(i,6)
                        if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)  %Z2->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)               %Z2->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP+b;
                            end
                        end
                    elseif XC >= Fibre_pos(i,6) && XC <= a-Fibre_pos(i,6) && YC > b-Fibre_pos(i,6)
                        if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) %Z3->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)            %Z3->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP-b;
                            end
                        end
                    elseif YC >= Fibre_pos(i,6) && YC <= b-Fibre_pos(i,6) && XC < Fibre_pos(i,6)
                        if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)     %Z4->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)            %Z4->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP+a; Fibre_pos(i+ics,3) = Y_TMP;
                            end
                        end
                    elseif YC >= Fibre_pos(i,6) && YC <= b-Fibre_pos(i,6) && XC > a-Fibre_pos(i,6)
                        if X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)     %Z5->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP > a-Fibre_pos(i,6)                %Z5->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP-a; Fibre_pos(i+ics,3) = Y_TMP;
                           end
                        end
                    end
                elseif Fibre_pos(i,4) == 4
                    if XC < Fibre_pos(i,6) && YC < Fibre_pos(i,6)
                        if X_TMP < Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)          %Z6->Z6
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP+a,Y_TMP+b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                        Fibre_pos(i+2,2) = X_TMP+a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP+a; Fibre_pos(i+3,3) = Y_TMP+b;
                                    end
                                end
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z6->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)             %Z6->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b; Fibre_pos(i+1,4) = 2;
                                N_fibre = N_fibre - 2; Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)            %Z6->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP+a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        end
                    elseif XC > a-Fibre_pos(i,6) && YC < Fibre_pos(i,6)
                        if X_TMP > a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)          %Z7->Z7
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP-a,Y_TMP+b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                        Fibre_pos(i+2,2) = X_TMP-a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP-a; Fibre_pos(i+3,3) = Y_TMP+b;
                                    end
                                end
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z7->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)             %Z7->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        elseif X_TMP > a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)            %Z7->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP-a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        end
                    elseif XC > a-Fibre_pos(i,6) && YC > b-Fibre_pos(i,6)
                        if X_TMP > a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)        %Z8->Z8
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP-a,Y_TMP-b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                        Fibre_pos(i+2,2) = X_TMP-a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP-a; Fibre_pos(i+3,3) = Y_TMP-b;
                                    end
                                end
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z8->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)             %Z8->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        elseif X_TMP > a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)            %Z8->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP-a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        end
                    elseif XC < Fibre_pos(i,6) && YC > b-Fibre_pos(i,6)
                        if X_TMP < Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)            %Z9->Z9
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP+a,Y_TMP-b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                        Fibre_pos(i+2,2) = X_TMP+a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP+a; Fibre_pos(i+3,3) = Y_TMP-b;
                                    end
                                end
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z9->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)             %Z9->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)              %Z9->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP+a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                            end
                        end
                    end
                end
            end
            if i>N_fibre, break,end
            if Fibre_pos(i,4) == 4
                fibre_jump = 3;
            end
        end
    end
    disp('Time elapsed after step 2 [min]: '); disp(toc/60);
    if image_option == 1
        index1 = 2;
        f_image_per(index1);
    end
    %
    % Second Heuristic - Compact the fibres on the outskirts of the RVE
    %
    GO2 = 0;
    if Square_size > a/2-Rmax, GO2 = 1; end
    if N_fibre <= Max_fibres && (mod(N_cycles,2) == 0 || GO2 == 1)
        fibre_jump = 0;
        for i=1:2:N_fibre
            if i>N_fibre, break, end
            if fibre_jump ~= 0
                fibre_jump = fibre_jump - 1;
                continue
            end
            % if fibre is on the LEFT side of the RVE... (Position 1)
            if Fibre_pos(i, 2) < Square_size && Fibre_pos(i,3) > Square_size &&...
                    Fibre_pos(i,3) < b-Square_size
                MIN = a*2; THETA_MIN = 3; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=-pi/2:pi/90:pi/2
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 3
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    else
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if X_TMP >= Fibre_pos(i,6)                       %Z4->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                            continue
                        else                                %Z4->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP+a; Fibre_pos(i+ics,3) = Y_TMP;
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the RIGHT side of the RVE... (Position 2)
            if Fibre_pos(i,2) > a-Square_size && Fibre_pos(i,3) > Square_size &&...
                    Fibre_pos(i,3) < b-Square_size
                MIN = a*2; THETA_MIN = 0; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=pi/2:pi/90:1.5*pi
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 0
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    else
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if X_TMP <= a-Fibre_pos(i,6)                     %Z5->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                            continue
                        else                                %Z5->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP-a; Fibre_pos(i+ics,3) = Y_TMP;
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the BOTTOM side of the RVE... (Position 3)
            if Fibre_pos(i,3) < Square_size && Fibre_pos(i,2) > Square_size &&...
                    Fibre_pos(i,2) < a-Square_size
                MIN = a*2; THETA_MIN = -1; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=0:pi/90:pi
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= -1
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    else
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Y_TMP >= Fibre_pos(i,6)                       %Z2->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                            continue
                        else                                %Z2->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP+b;
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the TOP side of the RVE... (Position 4)
            if Fibre_pos(i,3) > b-Square_size && Fibre_pos(i,2) > Square_size &&...
                    Fibre_pos(i,2) < a-Square_size
                MIN = a*2; THETA_MIN = 2; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=-pi:pi/90:0
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 2
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    else
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Y_TMP <= b-Fibre_pos(i,6)                     %Z3->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                            Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                            Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                            continue
                        else                                %Z3->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP-b;
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the LEFT-BOTTOM side of the RVE... (Position 5)
            if Fibre_pos(i,2) < Square_size && Fibre_pos(i,3) < Square_size
                MIN = a*2; THETA_MIN = 2; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=0:pi/90:pi/2
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 2
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    elseif Fibre_pos (i,4) == 2
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Fibre_pos(i,2) >= Fibre_pos(i,6)              %Z2
                            if Y_TMP >= Fibre_pos(i,6)                   %Z2->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                            %Z2->Z2
                                [check dist] = f_overlap(X_TMP,Y_TMP+b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP+b;
                                    continue
                                end
                            end
                        else                                 %Z4
                            if X_TMP >= Fibre_pos(i,6)                    %Z4->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                             %Z4->Z4
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP+a; Fibre_pos(i+ics,3) = Y_TMP;
                                    continue
                                end
                            end
                        end
                    else
                        if X_TMP < Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)            %Z6->Z6
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP+a,Y_TMP+b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                        Fibre_pos(i+2,2) = X_TMP+a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP+a; Fibre_pos(i+3,3) = Y_TMP+b;
                                        fibre_jump = 3; continue
                                    else fibre_jump = 3; continue
                                    end
                                else fibre_jump = 3; continue
                                end
                            else fibre_jump = 3; continue
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z6->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];  %% POSSIBLE BUG
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                            continue
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)             %Z6->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)            %Z6->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP+a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the RIGHT-BOTTOM side of the RVE... (Position 6)
            if Fibre_pos(i,2) > a-Square_size && Fibre_pos(i,3) < Square_size
                MIN = a*2; THETA_MIN = 0; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=pi/2:pi/90:pi
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 0
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    elseif Fibre_pos (i,4) == 2
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Fibre_pos(i,2) <= a-Fibre_pos(i,6)            %Z2
                            if Y_TMP >= Fibre_pos(i,6)                   %Z2->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                            %Z2->Z2
                                [check dist] = f_overlap(X_TMP,Y_TMP+b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP+b;
                                    continue
                                end
                            end
                        else                                 %Z5
                            if X_TMP <= a-Fibre_pos(i,6)                  %Z5->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                             %Z5->Z5
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP-a; Fibre_pos(i+ics,3) = Y_TMP;
                                    continue
                                end
                            end
                        end
                    else                                     %Z7
                        if X_TMP > a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)          %Z7->Z7
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP-a,Y_TMP+b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                        Fibre_pos(i+2,2) = X_TMP-a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP-a; Fibre_pos(i+3,3) = Y_TMP+b;
                                        fibre_jump = 3; continue
                                    else fibre_jump = 3; continue
                                    end
                                else fibre_jump = 3; continue
                                end
                            else fibre_jump = 3; continue
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z7->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                            continue
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)             %Z7->Z2
                            [check dist] = f_overlap(X_TMP,Y_TMP+b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP+b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        elseif X_TMP > a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)            %Z7->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP-a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the RIGHT-TOP side of the RVE... (Position 7)
            if Fibre_pos(i,2) > a-Square_size && Fibre_pos(i,3) > b-Square_size
                MIN = a*2; THETA_MIN = 0; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=-pi:pi/90:-pi/2
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*Fibre_pos(i,6); end
                end
                if THETA_MIN ~= 0
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    elseif Fibre_pos(i,4) == 2
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Fibre_pos(i,2) <= a-Fibre_pos(i,6)            %Z3
                            if Y_TMP <= b-Fibre_pos(i,6)                 %Z3->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                            %Z3->Z3
                                [check dist] = f_overlap(X_TMP,Y_TMP-b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP-b;
                                    continue
                                end
                            end
                        else                                 %Z5
                            if X_TMP <= a-Fibre_pos(i,6)                  %Z5->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                             %Z5->Z5
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP-a; Fibre_pos(i+ics,3) = Y_TMP;
                                    continue
                                end
                            end
                        end
                    else                                     %Z8
                        if X_TMP > a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)        %Z8->Z8
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP-a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP-a,Y_TMP-b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                        Fibre_pos(i+2,2) = X_TMP-a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP-a; Fibre_pos(i+3,3) = Y_TMP-b;
                                        fibre_jump = 3; continue
                                    else fibre_jump = 3; continue
                                    end
                                else fibre_jump = 3; continue
                                end
                            else fibre_jump = 3; continue
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z8->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                            continue
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)             %Z8->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        elseif X_TMP > a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)            %Z8->Z5
                            [check dist] = f_overlap(X_TMP-a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP-a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        end
                    end
                end
            end
            % if fibre is on the LEFT-TOP side of the RVE... (Position 8)
            if Fibre_pos(i,2) < Square_size && Fibre_pos(i,3) > b-Square_size
                MIN = a*2; THETA_MIN = 1; RAD = 0.75*Fibre_pos(i,6); GO = 0;
                while RAD ~= 0 && GO == 0
                    for j=-pi/2:pi/90:0
                        X_TMP = Fibre_pos(i,2) + RAD*cos(j);
                        Y_TMP = Fibre_pos(i,3) + RAD*sin(j);
                        [check dist] = f_overlap(X_TMP,Y_TMP,i,Fibre_pos(i,6),Fibre_pos(i,5));
                        if check == 0 && dist < MIN
                            GO = 1;
                            MIN = dist;
                            THETA_MIN = j;
                        end
                    end
                    if GO == 0, RAD = RAD - 0.25*R; end
                end
                if THETA_MIN ~= 1
                    X_TMP = Fibre_pos(i,2) + RAD*cos(THETA_MIN);
                    Y_TMP = Fibre_pos(i,3) + RAD*sin(THETA_MIN);
                    if Fibre_pos(i,4) == 0                  %Z1->Z1
                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                        continue
                    elseif Fibre_pos(i,4) == 2
                        if i ~= N_fibre
                            if Fibre_pos(i,2) == Fibre_pos(i+1,2) || Fibre_pos(i,3) == Fibre_pos(i+1,3)
                                ics = 1;
                            else
                                ics = -1;
                            end
                        else ics = -1;
                        end
                        if Fibre_pos(i,2) >= Fibre_pos(i,6)              %Z3
                            if Y_TMP <= b-Fibre_pos(i,6)                 %Z3->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                            %Z3->Z3
                                [check dist] = f_overlap(X_TMP,Y_TMP-b,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP; Fibre_pos(i+ics,3) = Y_TMP-b;
                                    continue
                                end
                            end
                        else                                    %Z4
                            if X_TMP >= Fibre_pos(i,6)                       %Z4->Z1
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                Fibre_pos(i,4) = 0; N_fibre = N_fibre - 1;
                                Fibre_pos(i+ics,:) = []; Vec_mem(i+ics,:) = [];
                                continue
                            else                                 %Z4->Z4
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+ics,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                    Fibre_pos(i+ics,2) = X_TMP+a; Fibre_pos(i+ics,3) = Y_TMP;
                                    continue
                                end
                            end
                        end
                    else
                        if X_TMP < Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)            %Z9->Z9
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                [check dist] = f_overlap(X_TMP+a,Y_TMP,i+2,Fibre_pos(i,6),Fibre_pos(i,5));
                                if check == 0
                                    [check dist] = f_overlap(X_TMP+a,Y_TMP-b,i+3,Fibre_pos(i,6),Fibre_pos(i,5));
                                    if check == 0
                                        Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP;
                                        Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                        Fibre_pos(i+2,2) = X_TMP+a; Fibre_pos(i+2,3) = Y_TMP;
                                        Fibre_pos(i+3,2) = X_TMP+a; Fibre_pos(i+3,3) = Y_TMP-b;
                                        fibre_jump = 3; continue
                                    else fibre_jump = 3; continue
                                    end
                                else fibre_jump = 3; continue
                                end
                            else fibre_jump = 3; continue
                            end
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6)    %Z9->Z1
                            Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 0;
                            N_fibre = N_fibre - 3;
                            Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = []; Fibre_pos(i+1,:) = [];
                            Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = []; Vec_mem(i+1,:) = [];
                            continue
                        elseif X_TMP >= Fibre_pos(i,6) && X_TMP <= a-Fibre_pos(i,6) && Y_TMP > b-Fibre_pos(i,6)             %Z9->Z3
                            [check dist] = f_overlap(X_TMP,Y_TMP-b,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP; Fibre_pos(i+1,3) = Y_TMP-b;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        elseif Y_TMP >= Fibre_pos(i,6) && Y_TMP <= b-Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)              %Z9->Z4
                            [check dist] = f_overlap(X_TMP+a,Y_TMP,i+1,Fibre_pos(i,6),Fibre_pos(i,5));
                            if check == 0
                                Fibre_pos(i,2) = X_TMP; Fibre_pos(i,3) = Y_TMP; Fibre_pos(i,4) = 2;
                                Fibre_pos(i+1,2) = X_TMP+a; Fibre_pos(i+1,3) = Y_TMP;
                                Fibre_pos(i+1,4) = 2; N_fibre = N_fibre - 2;
                                Fibre_pos(i+2,:) = []; Fibre_pos(i+2,:) = [];
                                Vec_mem(i+2,:) = []; Vec_mem(i+2,:) = [];
                                continue
                            end
                        end
                    end
                end
            end
            if Fibre_pos(i,4) == 4
                fibre_jump = 3;
            end
        end
        r=max(R,R1);
        if Square_size <= a/2 - r
            Square_size = Square_size + Square_inc;
        end
    end
    %
    % Are there any rogue fibres??? --> WASTE THEM!!
    %
    % Waste the random fibres that appear in the outskirts of the RVE 
    %    (POSSIBLE BUG in the previous routines)
    for i=1:1:N_fibre
        if i>=N_fibre, break, end
        if Fibre_pos(i,4)==0
           
           if Fibre_pos(i,2)<-Fibre_pos(i,6) ||Fibre_pos(i,2)>a+Fibre_pos(i,6)||...
              Fibre_pos(i,3)<-Fibre_pos(i,6) || Fibre_pos(i,3)>b+Fibre_pos(i,6)
          
              if Fibre_pos(i,5)==0
                     Vol_fibre = Vol_fibre - A_1_fibre/A_total;
                     Vol_fibre_1=Vol_fibre_1-A_1_fibre/A_total;
              else
                     Vol_fibre = Vol_fibre - A_2_fibre/A_total;
                     Vol_fibre_2=Vol_fibre_2-A_2_fibre/A_total;
              end
              Fibre_pos(i,:) = []; Vec_mem(i,:) = [];         
              N_fibre = N_fibre - 1;
              N_fibre_real = N_fibre_real - 1;
           end
        end
    
    end
    GO = 0; Lim = 0.90;
    while GO==0
        GO = 1;
        for i=1:1:N_fibre
            if i>=N_fibre, break, end
            if Fibre_pos(i,4) == 2
                X_TMP = Fibre_pos(i,2); Y_TMP = Fibre_pos(i,3);
                if X_TMP < -Lim*Fibre_pos(i,6) || (X_TMP > Lim*Fibre_pos(i,6) && X_TMP < Fibre_pos(i,6)) || ...
                        Y_TMP < -Lim*Fibre_pos(i,6) || (Y_TMP > Lim*Fibre_pos(i,6) && Y_TMP < Fibre_pos(i,6)) || ...
                        X_TMP > a+Lim*Fibre_pos(i,6) || (X_TMP > a-Fibre_pos(i,6) && X_TMP < a-Lim*Fibre_pos(i,6)) || ...
                        Y_TMP > b+Lim*Fibre_pos(i,6) || (Y_TMP > b-Fibre_pos(i,6) && Y_TMP < b-Lim*Fibre_pos(i,6))
                    if Fibre_pos(i,5)==0
                        Vol_fibre = Vol_fibre - A_1_fibre/A_total;
                        Vol_fibre_1=Vol_fibre_1-A_1_fibre/A_total;
                    else
                         Vol_fibre = Vol_fibre - A_2_fibre/A_total;
                        Vol_fibre_2=Vol_fibre_2-A_2_fibre/A_total;
                    end
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    N_fibre = N_fibre - 2;
                    N_fibre_real = N_fibre_real - 1;
                    GO = 0;
                    break
                end
            end
            if Fibre_pos(i,4) == 4
                X_TMP = Fibre_pos(i,2); Y_TMP = Fibre_pos(i,3);
                if (X_TMP>Lim*Fibre_pos(i,6) && X_TMP<Fibre_pos(i,6) && Y_TMP<Fibre_pos(i,6)) || ...
                        (Y_TMP>Lim*Fibre_pos(i,6) && Y_TMP<Fibre_pos(i,6) && X_TMP<Fibre_pos(i,6)) || ...
                        (X_TMP>a-Fibre_pos(i,6) && X_TMP<a-Lim*Fibre_pos(i,6) && Y_TMP<Fibre_pos(i,6)) || ...
                        (Y_TMP<Fibre_pos(i,6) && Y_TMP>Lim*Fibre_pos(i,6) && X_TMP>a-Fibre_pos(i,6)) || ...
                        (X_TMP>a-Fibre_pos(i,6) && X_TMP<a-Lim*Fibre_pos(i,6) && Y_TMP>b-Fibre_pos(i,6)) || ...
                        (Y_TMP>b-Fibre_pos(i,6) && Y_TMP<b-Lim*Fibre_pos(i,6) && X_TMP>a-Fibre_pos(i,6)) || ...
                        (X_TMP>Lim*Fibre_pos(i,6) && X_TMP<Fibre_pos(i,6) && Y_TMP>b-Fibre_pos(i,6)) || ...
                        (Y_TMP>b-Fibre_pos(i,6) && Y_TMP<b-Lim*Fibre_pos(i,6) && X_TMP<Fibre_pos(i,6))
                     if Fibre_pos(i,5)==0
                        Vol_fibre = Vol_fibre - A_1_fibre/A_total;
                        Vol_fibre_1=Vol_fibre_1-A_1_fibre/A_total;
                    else
                         Vol_fibre = Vol_fibre - A_2_fibre/A_total;
                        Vol_fibre_2=Vol_fibre_2-A_2_fibre/A_total;
                    end
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];                   %% POSSIBLE BUG
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    N_fibre = N_fibre - 4;
                    N_fibre_real = N_fibre_real - 1;
                   
                    GO = 0;
                    break
                end
                dist_tmp1 = sqrt(X_TMP^2+Y_TMP^2);
                dist_tmp2 = sqrt((a-X_TMP)^2+Y_TMP^2);
                dist_tmp3 = sqrt((a-X_TMP)^2+(b-Y_TMP)^2);
                dist_tmp4 = sqrt(X_TMP^2+(b-Y_TMP)^2);
                if (dist_tmp1 > Fibre_pos(i,6)-S_base/2 && dist_tmp1 < Fibre_pos(i,6)+S_base/2) || ...
                        (dist_tmp2 > Fibre_pos(i,6)-S_base/2 && dist_tmp2 < Fibre_pos(i,6)+S_base/2) || ...
                        (dist_tmp3 > Fibre_pos(i,6)-S_base/2 && dist_tmp3 < Fibre_pos(i,6)+S_base/2) || ...
                        (dist_tmp4 > Fibre_pos(i,6)-S_base/2 && dist_tmp4 < Fibre_pos(i,6)+S_base/2)
                     if Fibre_pos(i,5)==0
                        Vol_fibre = Vol_fibre - A_1_fibre/A_total;
                        Vol_fibre_1=Vol_fibre_1-A_1_fibre/A_total;
                    else
                         Vol_fibre = Vol_fibre - A_2_fibre/A_total;
                        Vol_fibre_2=Vol_fibre_2-A_2_fibre/A_total;
                    end
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];                   %% POSSIBLE BUG
                    Fibre_pos(i,:) = []; Vec_mem(i,:) = [];
                    N_fibre = N_fibre - 4;
                    N_fibre_real = N_fibre_real - 1;
                   
                    GO = 0;
                    break
                end
            end
        end
    end
    disp('Time elapsed after step 3 [min]: '); disp(toc/60);
    disp('Achieved Fibre Volume: '); disp(Vol_fibre);
    if image_option == 1
        index1 = 3;
        f_image_per(index1);
    end
    %
    N_cycles = N_cycles + 1;
    if N_cycles == N_cycles_max + 1, break, end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Miscellaneous Information                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Dif_error = 100*(Vol_fibre-Vol_fibre_req)/Vol_fibre_req;
if Dif_error < 0
    disp(' ');
    disp('Fibre Volume Fraction NOT REACHED');
    disp('Error in Volume Fibre [%]: '); disp(Dif_error);
    status = 0;
else
    disp(' ');
    disp('Random Distribution of Fibres COMPLETED');
    disp('Number of Attempts: '); disp(N_attempts);
    disp('Number of fibres: '); disp(N_fibre_real);
    disp('Elapsed Time [min]: '); disp(toc/60);
    disp('Achieved Fibre Volume: '); disp(Vol_fibre);
    disp('Achieved Type 1 Fibre Volume: '); disp(Vol_fibre_1);
    disp('Achieved Type 2 Fibre Volume: '); disp(Vol_fibre_2);
    disp('Error in Volume Fibre [%]: '); disp(Dif_error);
    status = 1;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
