function force_matrix = Acini_ForceGen_v2_T1(ECM, point_set, SD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization

force_matrix = zeros(200,1);   %Zeroes All Force Vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Force Generation

%for time_interval = 1:8000                                       %Time Interval For Force Generation
    
    %Cell 1
    %For Points = 1:1:20
    
    %point_set_1 = [1:20];                       %Points Able to Generate Force
    point_1 = randi(length(point_set{1}));                        %Random Place Number Generated
    force_location = point_set{1}(point_1);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{1},SD{1});    %Force Generation Value
    
    %Cell 2
    %For Points = 21:1:40
    
    %point_set_2 = [21:40];                  %Points Able to Generate Force
    point_2 = randi(length(point_set{2}));                        %Random Place Number Generated
    force_location = point_set{2}(point_2);                       %Random Point Selected From Set    
    force_matrix(force_location) = normrnd(ECM{2},SD{2});    %Force Generation Value
    
    %Cell 3
    %For Points = 41:1:60
     
    %point_set_3 = [41:60];                     %Points Able to Generate Force
    point_3 = randi(length(point_set{3}));                        %Random Place Number Generated
    force_location = point_set{3}(point_3);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{3},SD{3});    %Force Generation Value
    
    %Cell 4
    %For Points = 61:1:80
    
    %point_set_4 = [61:80];                        %Points Able to Generate Force
    point_4 = randi(length(point_set{4}));                        %Random Place Number Generated
    force_location = point_set{4}(point_4);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{4},SD{4});    %Force Generation Value
    
    %Cell 5
    %For Points = 81:1:100
    
    %point_set_5 = [81:100];                  %Points Able to Generate Force
    point_5 = randi(length(point_set{5}));                        %Random Place Number Generated
    force_location = point_set{5}(point_5);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{5},SD{5});    %Force Generation Value
    
    %Cell 6
    %For Points = 101:1:120
    
    %point_set_6 = [101:120];         %Points Able to Generate Force
    point_6 = randi(length(point_set{6}));                        %Random Place Number Generated
    force_location = point_set{6}(point_6);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{6},SD{6});    %Force Generation Value
    
    %Cell 7
    %For Points = 121:1:140
    
    %point_set_7 = [121:140];         %Points Able to Generate Force
    point_7 = randi(length(point_set{7}));                        %Random Place Number Generated
    force_location = point_set{7}(point_7);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{7},SD{7});    %Force Generation Value
    
    %Cell 8
    %For Points = 141:1:160
     
    %point_set_8 = [141:160];                 %Points Able to Generate Force
    point_8 = randi(length(point_set{8}));                        %Random Place Number Generated
    force_location = point_set{8}(point_8);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{8},SD{8});    %Force Generation Value

    %Cell 9
    %For Points = 161:1:180
    
    %point_set_9 = [161:180]; %Points Able to Generate Force
    point_9 = randi(length(point_set{9}));                        %Random Place Number Generated
    force_location = point_set{9}(point_9);                       %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{9},SD{9});    %Force Generation Value

    %Cell 10
    %For Points = 181:1:200
     
    %point_set_10 = [181:200];            %Points Able to Generate Force
    point_10 = randi(length(point_set{10}));                      %Random Place Number Generated
    force_location = point_set{10}(point_10);                     %Random Point Selected From Set
    force_matrix(force_location) = normrnd(ECM{10},SD{10});   %Force Generation Value
    
%end
end
