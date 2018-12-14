function dydt = Acini_Function_v1_T1(t,i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global Parameters

global ncell                %Number of Cells
global npoints              %Number of Cell Points
global c                    %Cell Number
global Cell_Point_Number    %Point on Cell

global l0_out               %Initial Length Between Cell Points, um
global l0_in                %Initial Length Between Nucleus Points, um
global l0_inout             %Initial Length Between Cell Points and Nucleus Points, um
global l0_nu                %Initial Length Between Internal Nucleus Points, um
global l_bond               %Initial Length Between Cell to Cell Binding, um
global thresh_max           %Max Cell to Cell binding Distance, um
global thresh_min           %Min Cell to Cell binding Distance, um
global dist_nr

global k_cout               %Stiffness Value of Cell Membrane, nN/um
global k_cin                %Stiffness Value of Nucleus Membrane, nN/um
global k_nu                 %Stiffness Value of Internal Nucleus, nN/um
global k_matrix             %Stiffness Value of ECM
global k_cell               %Stiffness Value of Cell to Cell Binding, nN/um
global k_in                 %Stiffness Value of Cell Membrane to Nucleus, nN/um

global F_external           %Forces on Cell Points, nN
global F_intx               %Forces on Nucleus Points X-direction, nN
global F_inty               %Forces on Nucleus Points X-direction, nN
global force_number

global v                    %Viscosity Coefficient, nN*s/um

global total_t              %Total Time of Simulation

global ecm_thresh_max       %Max Cell to ECM binding Distance, um
global ecm_thresh_min       %Max Cell to ECM binding Distance, um
global ecm_point            %ECM Point Number
global ecm_x                %ECM X-Position
global ecm_y                %ECM Y-Position
global ecm_bond             %Initial Length Between Cell to ECM Binding, um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization

attached =zeros(1,4*ncell*npoints);        %Zeros all Cell to Cell Adhesions Attachment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Force Calculation on Cell Points From Neighboring Cell Points and Nucleus Points

for s=1:(ncell*npoints)                %Iteration, Points 1 Through (ncell*npoints) of Cell Points
    
    
    if Cell_Point_Number(s)>1 && Cell_Point_Number(s)<20                                                                             %If Cell Point # 1<CP<npoints use this for loop
        d1(s)=sqrt((i(s)-i(s+1))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((i(s)-i(s-1))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
        F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(i(s+1)-i(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(i(s-1)-i(s))/d2(s);                       %Force Calculation Between Cell Points X-Direction
        F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(i(s+(ncell*npoints)+1)-i(s+(ncell*npoints)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
            k_cout(s)*(d2(s)-l0_out)*(i(s+(ncell*npoints)-1)-i(s+(ncell*npoints)))/d2(s);
        
    elseif Cell_Point_Number(s) == 1                                                                                                 %If Cell Point # = 1 use this for loop
        d1(s)=sqrt((i(s)-i(s+1))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((i(s)-i(s+(npoints-1)))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)+(npoints-1)))^2);                             %Distance from Cell Point to Cell Point Right
        F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(i(s+1)-i(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(i(s+(npoints-1))-i(s))/d2(s);             %Force Calculation Between Cell Points X-Direction
        F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(i(s+(ncell*npoints)+1)-i(s+(ncell*npoints)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
            k_cout(s)*(d2(s)-l0_out)*(i(s+(ncell*npoints)+(npoints-1))-i(s+(ncell*npoints)))/d2(s);
        
    elseif Cell_Point_Number(s) == npoints                                                                                           %If Cell Point # = npoints use this for loop
        d1(s)=sqrt((i(s)-i(s-(npoints-1)))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)-(npoints-1)))^2);                             %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((i(s)-i(s-1))^2+(i(s+(ncell*npoints))-i(s+(ncell*npoints)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
        F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(i(s-(npoints-1))-i(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(i(s-1)-i(s))/d2(s);
        F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(i(s+(ncell*npoints)-(npoints-1))-i(s+(ncell*npoints)))/d1(s)+...                       %Force Calculation Between Cell Points Y-Direction
            k_cout(s)*(d2(s)-l0_out)*(i(s+(ncell*npoints)-1)-i(s+(ncell*npoints)))/d2(s);
        
    end
    
    d0(s)=sqrt((i(s)-i(s+(2*ncell*npoints)))^2+(i(s+(ncell*npoints))-i(s+(3*ncell*npoints)))^2);       %Distance From Cell Point to Corresponding Nucleus Point
    
    F_outinx(s)=k_in(s)*(d0(s)-l0_inout)*(i(s+(2*ncell*npoints))-i(s))/d0(s);                          %Force Calculation Between Cell Point and Nucleus Point X-Direction
    F_outiny(s)=k_in(s)*(d0(s)-l0_inout)*(i(s+(3*ncell*npoints))-i(s+(ncell*npoints)))/d0(s);          %Force Calculation Between Cell Point and Nucleus POint Y-Direction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Initialization
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Force Calculation on Cell to Cell Adhesion
    
    F_adhx=0;                                   %Zeros all Cell to Cell Adhesions Forces
    F_adhy=0;                                   %Zeros all Cell to Cell Adhesions Forces
    
    if attached(s)==0                                                       %If Cell Point Is Not Attached Go Through This Loop
        for w=1:ncell*npoints                                               %Iteration, Points 1 through (ncell*npoints) of Cell Structure
            if c(w)~=c(s) && attached(w) == 0 && attached(s)==0             %If Points Not From Same Cell and Are Not Attached Proceed
                
                dist(s,w)=(sqrt((i(s)-i(w))^2+(i(s+ncell*npoints)-i(w+ncell*npoints))^2));                        %Distance From Cell Point To All Other Points of Sorrounding Cell Points
                if dist(s,w)<thresh_max % && dist(s,w)>thresh_min                                                   %If Cell Points Within Distance Threshold Proceed
                    attached(s)=w;                                                                                %Cell Point S Becomes Attached
                    attached(w)=s;                                                                                %Cell Point W Becomes Attached
                    F_adhx=F_adhx+k_cell*(dist(s,w)-l_bond)*(i(w)-i(s))/dist(s,w);                                %Adhesion Forces in X-direction
                    F_adhy=F_adhy+k_cell*(dist(s,w)-l_bond)*(i(w+ncell*npoints)-i(s+ncell*npoints))/dist(s,w);    %Adhesion Forces in Y-direction
                    
                else
                    F_adhx=F_adhx+0;                                        %Adhesion Force Unchanged
                    F_adhy=F_adhy+0;                                        %Adhesion Force Unchanged
                end
            else
                F_adhx=F_adhx+0;                                        %Adhesion Force Unchanged
                F_adhy=F_adhy+0;
            end
            %if c(w)~=c(s) && attached(w) == s && attached(s)==w
            %dist(s,w)=(sqrt((i(s)-i(w))^2+(i(s+ncell*npoints)-i(w+ncell*npoints))^2));
            %if dist(s,w)<new_cell_thresh && dist(s,w)>new_cell_thresh
            %add new cell
            %alter the npoints ncells in correct order
            %end
        end
        
    elseif attached(s) ~= 0
        w = attached(s);
        dist(s,w)=(sqrt((i(s)-i(w))^2+(i(s+ncell*npoints)-i(w+ncell*npoints))^2));                        %Distance From Cell Point To All Other Points of Sorrounding Cell Points
        if dist(s,w)<thresh_max  %&& dist(s,w)>thresh_min                                                   %If Cell Points Within Distance Threshold Proceed
            F_adhx=F_adhx+k_cell*(dist(s,w)-l_bond)*(i(w)-i(s))/dist(s,w);                                %Adhesion Forces in X-direction
            F_adhy=F_adhy+k_cell*(dist(s,w)-l_bond)*(i(w+ncell*npoints)-i(s+ncell*npoints))/dist(s,w);    %Adhesion Forces in Y-direction
            
        else
            attached(s)=0;
            attached(w)=0;
            F_adhx=F_adhx+0;                                        %Adhesion Force Unchanged
            F_adhy=F_adhy+0;                                        %Adhesion Force Unchanged
        end
        
        %if c(w)~=c(s) && attached(w) == s && attached(s)==w
        %dist(s,w)=(sqrt((i(s)-i(w))^2+(i(s+ncell*npoints)-i(w+ncell*npoints))^2));
        %if dist(s,w)<new_cell_thresh && dist(s,w)>new_cell_thresh
        %add new cell
        %alter the npoints ncells in correct order
        %end
        
        
    end
    
    F_adhx_out(s) = F_adhx;                                             %Set Adhesion Force
    F_adhy_out(s) = F_adhy;                                             %Set Adhesion Force
    
    if attached(s)~=0                                                   %If Attached Proceed
        F_adhx_out(attached(s)) = -F_adhx;                              %Opposite Reaction Force X-Direction
        F_adhy_out(attached(s)) = -F_adhy;                              %Opposite Reaction Force X-Direction
    end
    
    if attached(s) ~= 0
        
        F_extx (s) = 0;                                                 %Zero Force In X-Direction
        F_exty (s) = 0;                                                 %Zero Force In Y-Direction
        
    elseif attached(s) == 0
        
        
        r_x = mean(i(((2*ncell*npoints)+1):(3*ncell*npoints)));                              %Mean X-Position of Acini Structure
        r_y = mean(i(((3*ncell*npoints)+1):(4*ncell*npoints)));                              %Mean Y-Position of Acini Structure
        dist_center(s) = sqrt((i(s)-r_x)^2+(i(s+ncell*npoints)-r_y)^2);                      %Distance of Cell Point From Center of Acini Structure
        
        if c(s) == 1
            
            nr_x = mean(i(((2*ncell*npoints)+1):(2*ncell*npoints+npoints)));
            nr_y = mean(i(((3*ncell*npoints)+1):(3*ncell*npoints+npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 2
            
            nr_x = mean(i(((2*ncell*npoints)+npoints+1):(2*ncell*npoints+2*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+npoints+1):(3*ncell*npoints+2*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 3
            
            nr_x = mean(i(((2*ncell*npoints)+2*npoints+1):(2*ncell*npoints+3*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+2*npoints+1):(3*ncell*npoints+3*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 4
            
            nr_x = mean(i(((2*ncell*npoints)+3*npoints+1):(2*ncell*npoints+4*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+3*npoints+1):(3*ncell*npoints+4*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 5
            
            nr_x = mean(i(((2*ncell*npoints)+4*npoints+1):(2*ncell*npoints+5*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+4*npoints+1):(3*ncell*npoints+5*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 6
            
            nr_x = mean(i(((2*ncell*npoints)+5*npoints+1):(2*ncell*npoints+6*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+5*npoints+1):(3*ncell*npoints+6*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 7
            
            nr_x = mean(i(((2*ncell*npoints)+6*npoints+1):(2*ncell*npoints+7*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+6*npoints+1):(3*ncell*npoints+7*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 8
            
            nr_x = mean(i(((2*ncell*npoints)+7*npoints+1):(2*ncell*npoints+8*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+7*npoints+1):(3*ncell*npoints+8*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 9
            
            nr_x = mean(i(((2*ncell*npoints)+8*npoints+1):(2*ncell*npoints+9*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+8*npoints+1):(3*ncell*npoints+9*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        elseif c(s) == 10
            
            nr_x = mean(i(((2*ncell*npoints)+9*npoints+1):(2*ncell*npoints+10*npoints)));
            nr_y = mean(i(((3*ncell*npoints)+9*npoints+1):(3*ncell*npoints+10*npoints)));
            dist_n(s) = sqrt((i(s)-nr_x)^2+(i(s+ncell*npoints)-nr_y)^2);
            
        end
               
        if F_external(s) ~= 0  && dist_center(s) > dist_nr                                    %External Force Generated
            F_extx (s) = F_external(s)*(i(s)-nr_x)/dist_n(s);
            F_exty (s) = F_external(s)*(i(s+(ncell*npoints))-nr_y)/dist_n(s);
        else
            
            F_extx (s) = 0; %Zero Force In X-Direction
            F_exty (s) = 0; %Zero Force In Y-Direction
        end
    else
        
        F_extx (s) = 0; %Zero Force In X-Direction
        F_exty (s) = 0; %Zero Force In Y-Direction
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Total Cell Point Force Calculation
    
    Fx_out(s)= F_neighx(s)+F_outinx(s)+F_adhx_out(s)+F_extx(s);    %Sum Of Forces In X-Direction
    Fy_out(s)= F_neighy(s)+F_outiny(s)+F_adhy_out(s)+F_exty(s);    %Sum Of Forces in Y-Direction
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Force Calculation On Nucleus Points From Neighboring Nucleus Points and Cell Points

for q=1:(ncell*npoints) %Iteration, Points 1 Through (ncell*npoints) of Nucleus
    
    if  Cell_Point_Number(q) > 1 && Cell_Point_Number(q) < 20                                                               %If Nucleus Point # 1 < CP < npoints Use This For Loop
        
        d1_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)+1))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)+1))^2);       %Distance from Cell Point To Nucleus Point Left
        d2_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)-1))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)-1))^2);       %Distance from Cell Point To Nucleus Point Right
        F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)+1)-i(q+(2*ncell*npoints)))/d1_in(q)+...                                      %Force Calculation Between Nucleus Points X-Direction
            k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)-1)-i(q+(2*ncell*npoints)))/d2_in(q);
        F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)+1)-i(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)+...      %Force Calculation Between Nucleus Points Y-Direction
            k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)-1)-i(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
        
    elseif Cell_Point_Number(q) == 1                                                                                        %If Nucleus Point # = 1 Use This For Loop
        
        d1_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)+1))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)+1))^2);                     %Distance from Cell Point To Nucleus Point Left
        d2_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)+(npoints-1)))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)+(npoints-1)))^2); %Distance from Cell Point To Nucleus Point Right
        F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)+1)-i(q+(2*ncell*npoints)))/d1_in(q)+...                                                    %Force Calculation Between Nucleus Points X-Direction
            k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(npoints-1))-i(q+(2*ncell*npoints)))/d2_in(q);
        F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)+1)-i(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)+...                    %Force Calculation Between Nucleus Points Y-Direction
            k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)+(npoints-1))-i(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
        
    elseif Cell_Point_Number(q) == npoints                                                                                  %If Nucleus Point # = npoints Use This For Loop
        
        d1_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)-(npoints-1)))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)-(npoints-1)))^2); %Distance from Cell Point To Nucleus Point Left
        d2_in(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)-1))^2+(i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)-1))^2);                     %Distance from Cell Point To Nucleus Point Right
        F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)-(npoints-1))-i(q+(2*ncell*npoints)))/d1_in(q)+...                                          %Force Calculation Between Nucleus Points X-Direction
            k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)-1)-i(q+(2*ncell*npoints)))/d2_in(q);
        F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)-(npoints-1))-i(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)...           %Force Calculation Between Nucleus Points Y-Direction
            +k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(i(q+(2*ncell*npoints)+(ncell*npoints)-1)-i(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
        
    end
    
    d0(q)=sqrt((i(q)-i(q+(2*ncell*npoints)))^2+(i(q+(ncell*npoints))-i(q+(3*ncell*npoints)))^2);           %Distance From Cell Point to Corresponding Nucleus Point
    
    F_inoutx(q)=k_in(q)*(d0(q)-l0_inout)*(i(q)-i(q+(2*ncell*npoints)))/d0(q);                              %Force Calculation Between Cell Point and Nucleus Point X-Direction
    F_inouty(q)=k_in(q)*(d0(q)-l0_inout)*(i(q+(ncell*npoints))-i(q+(3*ncell*npoints)))/d0(q);              %Force Calculation Between Cell Point and Nucleus Point Y-Direction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Force Calculation On Internal Nucleus
    
    
    if  Cell_Point_Number(q) >= 1 && Cell_Point_Number(q) <= 10                                                         %If Nucleus Point Is >=1 & <=10 Use This For Loop
        
        d_nu(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)+10))^2+...                                           %Distance calculation Between Nucleus Points
            (i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)+10))^2);
        
        F_nx(q)=k_nu(q+(2*ncell*npoints))*(d_nu(q)-l0_nu)*(i(q+(2*ncell*npoints)+10)-i(q+(2*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
        F_ny(q)=k_nu(q+(3*ncell*npoints))*(d_nu(q)-l0_nu)*(i(q+(3*ncell*npoints)+10)-i(q+(3*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
        
    elseif Cell_Point_Number(q) >= 11 && Cell_Point_Number(q) <= 20                                                     %If Nucleus Point Is >=11 & <=20 Use This For Loop
        
        d_nu(q)=sqrt((i(q+(2*ncell*npoints))-i(q+(2*ncell*npoints)-10))^2+...                                           %Distance calculation Between Nucleus Points
            (i(q+(2*ncell*npoints)+(ncell*npoints))-i(q+(2*ncell*npoints)+(ncell*npoints)-10))^2);
        
        F_nx(q)=k_nu(q+(2*ncell*npoints))*(d_nu(q)-l0_nu)*(i(q+(2*ncell*npoints)-10)-i(q+(2*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
        F_ny(q)=k_nu(q+(3*ncell*npoints))*(d_nu(q)-l0_nu)*(i(q+(3*ncell*npoints)-10)-i(q+(3*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
        
    end
    
    F_inadhx=0;       %Zero Force In X-Direction
    F_inadhy=0;       %Zero Force In Y-Direction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Total Nucleus Force Calculation
    
    Fx_in(q)= F_inneighx(q)+F_inoutx(q)+F_inadhx+F_intx(q)+F_nx(q);    %Sum Of Forces In X-Direction
    Fy_in(q)= F_inneighy(q)+F_inouty(q)+F_inadhy+F_inty(q)+F_ny(q);    %Sum Of Forces In X-Direction
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization

dydt = zeros(ncell*npoints, 1);                     %Zero Integration Solutions

for d = 1:ncell*npoints                             %Iteration, Points 1 Through (ncell*npoints) of Cell Points
    
    dydt(d) = Fx_out(d)/v(d).';                     %Solves New Cell Point Position For Points 1 Through npoints in X-Direction
    dydt(d+ncell*npoints) = Fy_out(d)/v(d).';       %Solves New Cell Point Position For Points 1 Through npoints in Y-Direction
    dydt(d+2*ncell*npoints) = Fx_in(d)/v(d).';      %Solves New Nucleus Position For Points 1 Through npoints in X-Direction
    dydt(d+3*ncell*npoints) = Fy_in(d)/v(d).';      %Solves New Nucleus Position For Points 1 Through npoints in Y-Direction
    
end
end