clc
clear all

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

%Acini Structure Geometry

center = [0,0];  %Center of Structure
ncell = 10;      %Number of Nucleus on Structure
npoints = 20;    %Number of Points per Nucleus
r = 1;           %Radius of Nucleus, um
r1 = 17.5;       %Nucleus Array Radius, um

ncell_out = ncell;      %Number of Cells on Structure
npoints_out = npoints;  %Number of Points per Cell
r_out = 5;              %Radius of Cell, um
r2 = 17.5;              %Cell Array Radius, um

%Initialization for Acini Creation

point = 0;         %Starting Cell Point
point_out = 0;     %Starting Nucleus Point

for cell=1:ncell                                                                                       %For Loop of Cell Nucleus
    alpha=cell*2*pi/ncell;                                                                             %Angle of Cell Nucleus With Respect To Center of Structure
    cell_center(cell,:)=[r1*cos(alpha),r1*sin(alpha)];                                                 %Center of Cell Nucleus on Structure
    for boundary=1:npoints                                                                             %For loop of Points on Cell Nucleus
        point=point+1;                                                                                 %Selecting Next point
        theta=boundary*2*pi/npoints;                                                                   %Angle of Cell Nucleus Point With Respect To Center of Cell Nucleus
        pos(cell,boundary,:)=[cell_center(cell,1)+r*cos(theta), cell_center(cell,2)+r*sin(theta)];     %Positioning of Cell Nucleus
        cellpoint(point,:)= [pos(cell,boundary,1),pos(cell,boundary,2),cell,boundary];                 %Positioning of Cell Nucleus Points
    end
end

for cell_out=1:ncell_out                                                                                           %For Loop of Cells
    alpha_out=cell_out*2*pi/ncell_out;                                                                             %Angle of Cell With Respect To Center of Structure
    cell_center_out(cell_out,:)=[r2*cos(alpha_out),r2*sin(alpha_out)];                                             %Center of Cell on Structure
    for boundary_out=1:npoints_out                                                                                 %For loop of Points on Cell
        point_out=point_out+1;                                                                                     %Selecting Next point
        theta_out=boundary_out*2*pi/npoints_out;                                                                   %Angle of Cell Point With Respect To Center of Cell
        pos_out(cell_out,boundary_out,:) = ...                                                                     %Positioning of Cell
            [cell_center_out(cell_out,1)+r_out*cos(theta_out), cell_center_out(cell_out,2)+r_out*sin(theta_out)];
        cellpoint_out(point_out,:) = ...                                                                           %Positioning of Cell Points
            [pos_out(cell_out,boundary_out,1),pos_out(cell_out,boundary_out,2),cell_out,boundary_out];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ECM Structure

ecm_point = 100;                                                    %Number of ECM Points
ecm_r1 = 25; ecm_r2 = 50;                                           %Area of Point Location
ecm_r = sqrt(ecm_r1^2+(ecm_r2^2-ecm_r1^2)*rand(1,ecm_point));       %Using square root here ensures distribution uniformity by area
tau = 2*pi*rand(1,ecm_point);                                       %Angle of points on area
ecm_x = round(ecm_r.*cos(tau),0);                                   %ECM X-Position
ecm_y = round(ecm_r.*sin(tau),0);                                   %ECM Y-Position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Parameters

l0_out = sqrt((cellpoint_out(1,1)-cellpoint_out(2,1))^2+(cellpoint_out(1,2)-cellpoint_out(2,2))^2);  %Initial Length Between Cell Points, um
l0_in = sqrt((cellpoint(1,1)-cellpoint(2,1))^2+(cellpoint(1,2)-cellpoint(2,2))^2);                   %Initial Length Between Nucleus Points, um
l0_inout = sqrt((cellpoint_out(1,1)-cellpoint(1,1))^2+(cellpoint_out(1,2)-cellpoint(1,2))^2);        %Initial Length Between Cell Points and Nucleus Points, um

l0_nu = 4;                                       %Initial Length Between Internal Nucleus Points, um
l_bond = 0.25;                                   %Initial Length Between Cell to Cell Binding, um
ecm_bond = 0;                                    %Initial Length Between Cell to ECM Binding, um

dist_nr=22;

k_cout = 10 * ones(1,4*ncell*npoints);           %Stiffness Value of Cell Membrane, nN/um
k_cin  = 28  * ones(1,4*ncell*npoints);          %Stiffness Value of Nucleus Membrane, nN/um
k_in = 10 * ones(1,4*ncell*npoints);               %Stiffness Value of Cell Membrane to Nucleus, nN/um
k_nu = 28  * ones(1,4*ncell*npoints);            %Stiffness Value of Internal Nucleus, nN/um

est_k = 5;
k_in_new = 28 * ones(1,4*ncell*npoints);

k_matrix = 10;                                   %Stiffness Value of ECM
k_cell = 10;                                     %Stiffness Value of Cell to Cell Binding, nN/um

thresh_min = 0.25;                               %Min Cell to Cell binding Distance, um
thresh_max = 1;                                %Max Cell to Cell binding Distance, um
ecm_thresh_max = 7;                              %Max Cell to ECM binding Distance, um
ecm_thresh_min = 1;                              %Min Cell to ECM binding Distance, um


c = cellpoint_out(:,3);                                                           %Cell Number
Cell_Point_Number = ...
    vertcat(cellpoint_out(:,4),cellpoint_out(:,4),cellpoint(:,4),cellpoint(:,4)); %Point on Cell
v = 20 * ones(1,ncell*npoints);                                                   %Viscosity Coefficient, nN*s/um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulation Parameters
simulation_length = 30;  %Length of Simulation, s
time_step = 6;             %Time Step of Simulation, s
first_step = 6;            %Initial Time Step, s

%Initialization for Solver
Fx_in_matrix = zeros(ncell*npoints, simulation_length/time_step);    %Zeros Force Internal Matrix in X-direction, nN
Fy_in_matrix = zeros(ncell*npoints, simulation_length/time_step);    %Zeros Force Internal Matrix in y-direction, nN
Fx_out_matrix = zeros(ncell*npoints, simulation_length/time_step);   %Zeros Force External Matrix in X-direction, nN
Fy_out_matrix = zeros(ncell*npoints, simulation_length/time_step);   %Zeros Force External Matrix in Y-direction, nN
attached = zeros(1,4*ncell*npoints);                                 %Zeros all Cell to Cell Adhesions
F_adh_out = zeros(1,4*ncell*npoints);
F_intx = zeros(ncell*npoints, simulation_length/time_step);          %Zeros Force Internal in X-direction, nN
F_inty = zeros(ncell*npoints, simulation_length/time_step);          %Zeros Force Internal in Y-direction, nN
F_external = zeros(ncell*npoints,1);  %, simulation_length/time_step);      %Zeros Force External, nN
p = zeros(ncell*npoints*4, simulation_length/time_step*4);           %Zeros Point Positions, um
new_point = zeros(ncell*npoints*4, simulation_length/time_step*4);   %Zeros Point Positions, um
interval=0;                                                          %Initial Interval
force_number = 1;                                                    %Initial Force Vector
start = 0;                                                           %Start Time, s

%ECM and Force Generation On Selected Points
point_set_1 = [1:20];
point_set_2 = [21:40];
point_set_3 = [41:60];
point_set_4 = [61:80];
point_set_5 = [81:100];
point_set_6 = [101:120];
point_set_7 = [121:140];
point_set_8 = [141:160];
point_set_9 = [161:180];
point_set_10 = [181:200];
point_set = {point_set_1,point_set_2,point_set_3,point_set_4,point_set_5,point_set_6,point_set_7,point_set_8,point_set_9,point_set_10};

ECM1 = 10;                        %Random Force Generation Range for Cell 1, nN
SD1  = 2;
ECM2 = 10;                        %Random Force Generation Range for Cell 2, nN
SD2  = 2;
ECM3 = 10;                        %Random Force Generation Range for Cell 3, nN
SD3  = 2;
ECM4 = 10;                        %Random Force Generation Range for Cell 4, nN
SD4  = 2;
ECM5 = 10;                        %Random Force Generation Range for Cell 5, nN
SD5  = 2;
ECM6 = 10;                        %Random Force Generation Range for Cell 6, nN
SD6  = 2;
ECM7 = 10;                        %Random Force Generation Range for Cell 7, nN
SD7  = 2;
ECM8 = 10;                        %Random Force Generation Range for Cell 8, nN
SD8  = 2;
ECM9 = 10;                        %Random Force Generation Range for Cell 9, nN
SD9  = 2;
ECM10 = 10;                       %Random Force Generation Range for Cell 10, nN
SD10  = 2; 
ECM = {ECM1,ECM2,ECM3,ECM4,ECM5,ECM6,ECM7,ECM8,ECM9,ECM10};
SD  = {SD1, SD2, SD3 , SD4, SD5, SD6, SD7, SD8, SD9, SD10};

F_external  = Acini_ForceGen_v2_T1(ECM,point_set,SD); %Force External, nN

x0 = [19.5730684813128;18.7240969389947;17.3771928234621;15.6809873424956;13.8153334009029;11.9325065886132;10.5949350532538;9.50419651048900;8.53437794267337;7.91952470399752;8.40078328135734;9.30224124593738;10.6766390791327;12.3236256708423;14.0640866424636;15.7654765805010;17.4171456133828;18.8656966723753;19.5744549224811;19.8492024373430;11.0368897270638;10.1588737039425;8.83462860407123;7.15324199282598;5.28146030910109;3.41427065957708;1.71805151395631;0.272831801839068;-0.169538109725640;-0.239543257460237;-0.157932851575269;0.302910609655720;1.80205617471543;3.55305688674927;5.48488969130921;7.42964042692711;8.79013350538465;9.89210257798151;10.8696964761070;11.4917036148047;0.170026320629151;-0.272198848816417;-1.71715401944085;-3.41299074623717;-5.27992760366417;-7.15169262436917;-8.83331025545147;-10.1579543570738;-11.0364152100692;-11.4916692970801;-10.8698559738267;-9.89245111234668;-8.79054680812655;-7.42996824019207;-5.48525709954325;-3.55329885608455;-1.80200236170627;-0.302623357452272;0.158281485271702;0.239942899025338;-8.53470711090586;-9.50448076954549;-10.5950987121735;-11.9325915147584;-13.8154579028016;-15.6811472139658;-17.3774048616291;-18.7243621984486;-19.5733666440923;-19.8494987341219;-19.5747206662117;-18.8659597740174;-17.4174075892095;-15.7657624587588;-14.0643929953701;-12.3239306825332;-10.6769410099064;-9.30253582344239;-8.40108261225282;-7.91985157498942;-11.7901366629537;-12.5108746491018;-13.9666126132200;-15.6175632573223;-17.3124178401371;-19.0342062885685;-20.6355781465824;-21.9581415892703;-22.8355807363775;-23.1432231818540;-22.8355302109924;-21.9580525016546;-20.6354684191362;-19.0340879043944;-17.3122907353346;-15.6174321580408;-13.9664966325160;-12.5107837208746;-11.7900870717121;-11.4977516368974;-8.40105558835313;-9.30248221118050;-10.6768674365258;-12.3238398152569;-14.0642833026780;-15.7656428844382;-17.4172997902031;-18.8658655945455;-19.5746688535635;-19.8494683960804;-19.5733600925842;-18.7243877197037;-17.3774644517432;-15.6812325302633;-13.8155535229690;-11.9326836101389;-10.5951822552307;-9.50454376164211;-8.53474825097783;-7.91987235717267;0.158197685703508;-0.302698030046769;-1.80208472080177;-3.55337227097162;-5.48532148707334;-7.43002234782077;-8.79061324018733;-9.89253672246066;-10.8699516628012;-11.4917720291508;-11.0365196087705;-10.1580530242033;-8.83340656200311;-7.15178800876356;-5.28002181964115;-3.41308166719911;-1.71724154109438;-0.272286578138850;0.169943940005950;0.239858311214762;10.8697130072044;9.89213703801093;8.79017875399908;7.42968211975906;5.48492109883872;3.55307000756910;1.80204263403607;0.302862371629757;-0.157994608989334;-0.239608460308926;-0.169610126410688;0.272747820222531;1.71794999648085;3.41415144816114;5.28133659514255;7.15312815777350;8.83453365708824;10.1588006436129;11.0368421069306;11.4917022670856;19.5745073135344;18.8657454670879;17.4171896594399;15.7655194666459;14.0641274280782;12.3236653602358;10.6766817486959;9.30229098160102;8.40084083124015;7.91957981354415;8.53443326694879;9.50424430666873;10.5949699378329;11.9325373363816;13.8153832938815;15.6810528395400;17.3772693136048;18.7241770882445;19.5731439866666;19.8492668467616;22.8353340075154;21.9578768831917;20.6353063320498;19.0339202012151;17.3121196386118;15.6172473221607;13.9662949600119;12.5105715038005;11.7898661512673;11.4975214067089;11.7898863198884;12.5106090595617;13.9663424450401;15.6173013642103;17.3121761629574;19.0339760226534;20.6353587819098;21.9579182811650;22.8353570061257;23.1430088189347;11.9107272512940;13.6029685915076;14.9335273376936;15.7752739185713;16.0927725594843;15.9593583190936;14.8709531403025;13.4887576711384;12.0150236056592;10.3848028527422;8.49667137615919;6.75849712926098;5.36386751237904;4.32256547982180;4.40344132237891;4.87025197600829;5.48033269516710;6.41593458573034;8.16517504502049;10.0374698929580;18.1108070792054;19.7869310929111;21.1401933030157;22.0104153757989;22.2953771959588;21.9795141513788;21.1398538603150;19.9253688129422;18.2587904785488;16.4996203475222;14.7372424064231;13.0568117383294;11.8120311893192;10.9354596704966;10.6148475648690;10.7405131327820;11.8288890337013;13.2065608660140;14.6709439799765;16.2792872290064;18.2585048566440;19.9251146881022;21.1393285260618;21.9790904359773;22.2952807489804;22.0107456563911;21.1408967145498;19.7878298982425;18.1116828015402;16.2799781406364;14.6716931461206;13.2073052190867;11.8294594091501;10.7408928355265;10.6149454944797;10.9352510508790;11.8116802570778;13.0566336105440;14.7369971760393;16.4993406977286;12.0154218131954;13.4892161322013;14.8713933497687;15.9597766488837;16.0929663044647;15.7754504942884;14.9337120398814;13.6031475799633;11.9108808833526;10.0375902053190;8.16527011941149;6.41600054356987;5.48043157331431;4.87037841493987;4.40356379575041;4.32268461554632;5.36402390610455;6.75869751014993;8.49693306104106;10.3851401992996;1.93632137499141;3.74688064162866;4.70443006350507;5.32662286576622;5.80363900785328;5.89795909381208;4.89910439918202;3.54572313079235;1.86787048056697;-0.000129691003965956;-1.86811527398284;-3.54592747936925;-4.89924954949455;-5.89802462123512;-5.80367450973418;-5.32664095721693;-4.70444094125636;-3.74688507465463;-1.93630137907166;1.53218512283843e-05;-8.49696680390322;-6.75873623972108;-5.36405786720422;-4.32270016617383;-4.40358066233588;-4.87039030112901;-5.48044692885544;-6.41602323854326;-8.16525427298591;-10.0375521491113;-11.9108293503469;-13.6030962132803;-14.9336789097165;-15.7754512346091;-16.0930071834636;-15.9598347458714;-14.8714614673712;-13.4892787966011;-12.0154627811343;-10.3851665593923;-14.7369549641497;-13.0566039197797;-11.8116648075681;-10.9352300676104;-10.6149187702691;-10.7408909582186;-11.8294696905003;-13.2073349550308;-14.6717349288311;-16.2800169952451;-18.1117086891875;-19.7878378730274;-21.1408874925515;-22.0107217543743;-22.2952464752425;-21.9790492389268;-21.1392774507399;-19.9250505349590;-18.2584386415979;-16.4992810936759;-14.6709733003367;-13.2065930452885;-11.8289157629657;-10.7405330165050;-10.6148516691817;-10.9354463417148;-11.8120052255175;-13.0567728404692;-14.7371892100054;-16.4995504250333;-18.2587110121307;-19.9252911042699;-21.1397831512090;-21.9794707486537;-22.2953671814718;-22.0104335269487;-21.1402271578217;-19.7869680166150;-18.1108396229835;-16.2793139793050;-8.16513255647678;-6.41589567625772;-5.48030598543779;-4.87023618211832;-4.40343313007228;-4.32256491565273;-5.36387411904148;-6.75851260141813;-8.49669200051797;-10.3848231035282;-12.0150429492086;-13.4887742876635;-14.8709659889786;-15.9593681943994;-16.0927630052320;-15.7752588846722;-14.9335044553445;-13.6029359473421;-11.9106861572313;-10.0374251261874;1.86804228419993;3.54586726416382;4.89920532460965;5.89797086169256;5.80360460674712;5.32655557849384;4.70434724496703;3.74680373203873;1.93626683688168;-7.74859267661706e-06;-1.93627814273218;-3.74680521923120;-4.70434154251787;-5.32654019956004;-5.80357597845321;-5.89792547111241;-4.89913649364639;-3.54577778354140;-1.86793837275220;5.45680316998653e-05;15.8438071307067;15.5705836184806;15.1393754763921;14.6001137697998;14.0136988546975;13.4152395321248;12.9242708400749;12.5454634044385;12.2178062539734;12.0413698805557;12.1550597203994;12.4400348265489;12.8654844549622;13.3993030221387;14.0121578325616;14.6153132909056;15.1977496910421;15.6845373267063;15.9087545605020;15.9439368961202;7.24659219018414;6.91846687077973;6.48436658776204;5.94425683838718;5.34418274361472;4.74338933734947;4.17931542185338;3.72029369468248;3.53125607112781;3.51357346341294;3.55977422000094;3.76303773519223;4.20736617157362;4.75566631830681;5.35723367355324;5.96094454956135;6.46410403286169;6.86384341966872;7.21110833996474;7.39334105476972;-3.53219706859614;-3.72245506003511;-4.18245376700035;-4.74711705931019;-5.34818665175415;-5.94808996186102;-6.48756741168232;-6.92067704731003;-7.24759009583966;-7.39294835586573;-7.20932458876459;-6.86083755698721;-6.46015008548999;-5.95639078828061;-5.35245163962338;-4.75110073100737;-4.20345923299595;-3.76014020799981;-3.55806099723916;-3.51317753077218;-12.2181110017848;-12.5457602942673;-12.9245675497452;-13.4155555269799;-14.0140373278077;-14.6004580859846;-15.1397216067966;-15.5709254908134;-15.8441346319153;-15.9442364636373;-15.9090146424021;-15.6847612830080;-15.1979519284555;-14.6155144034203;-14.0123695865793;-13.3995241234426;-12.8657210894828;-12.4402947682958;-12.1553446597421;-12.0416733534414;-15.4624490781460;-15.7221129758504;-16.2099253799944;-16.7852367264216;-17.3794001072966;-17.9838081394104;-18.4923803746105;-18.8766730075926;-19.1481684619654;-19.2458717168771;-19.1482520299586;-18.8768283840380;-18.4925925305332;-17.9840556993781;-17.3796595629795;-16.7854793873994;-16.2101273649884;-15.7222565578552;-15.4625241343318;-15.3756532633927;-12.1553848732778;-12.4403915662488;-12.8658641471569;-13.3996968111895;-14.0125544712848;-14.6156920511590;-15.1981039459814;-15.6848685095360;-15.9090646080289;-15.9442142805809;-15.8440342765694;-15.5707554382314;-15.1395001455893;-14.6002070495768;-14.0137746425844;-13.4153056317611;-12.9243562118338;-12.5456077683993;-12.2180268997882;-12.0416546087673;-3.55815039942285;-3.76022008551277;-4.20353265299736;-4.75116746839884;-5.35251377587945;-5.95644960203775;-6.46020906971457;-6.86090315549168;-7.20940511545563;-7.39305127839638;-7.24771811158469;-6.92082611230308;-6.48773419110067;-5.94826645084541;-5.34836417999942;-4.74728630307195;-4.18260860183171;-3.72259189402268;-3.53231406691525;-3.51327888911248;7.21104430370370;6.86375041670681;6.46398670971865;5.96080889830259;5.35708894468853;4.75552387138911;4.20723591027832;3.76292701869997;3.55969175307922;3.51353135132402;3.53125543416991;3.72032954061145;4.17937738733396;4.74346561609256;5.34426456181803;5.94433239193938;6.48442246544092;6.91849230436905;7.24658616286911;7.39330637967701;15.9088232326191;15.6846114787397;15.1978287676570;14.6153967757736;14.0122439680131;13.3993885340633;12.8655671570775;12.4401130366545;12.1551313108903;12.0414324726237;12.2178609208255;12.5455116889815;12.9243148102299;13.4152830626717;14.0137444021797;14.6001617245172;15.1394266734084;15.5706393816524;15.8438677348313;15.9440012696559;19.1480120786921;18.8765717126792;18.4923218757790;17.9837717879997;17.3793713341766;16.7851981991571;16.2098595458896;15.7220029246489;15.4622850774377;15.3754361362142;15.4622607822249;15.7219568442276;16.2097952408206;16.7851213700053;17.3792892412232;17.9836934251586;18.4922553928262;18.8765235764533;19.1479863817071;19.2456572201120;10.7709081463308;11.3128793127166;11.7400074366476;12.0202198436900;12.1734674532174;12.1318770812491;11.7793545714507;11.2899023998921;10.7621279092917;10.1751103693452;9.57701323745120;9.03728766495208;8.59968646782232;8.30700993792951;8.29368848074287;8.44248424829870;8.65174867564124;9.00994485988499;9.56638268627292;10.1714835532148;17.0933625808180;17.6030542803638;18.0284302543594;18.3053597897076;18.3963886425429;18.3050522029632;18.0833547442817;17.6971559845154;17.1229060125118;16.5041549487749;15.8847643231847;15.3066404891902;14.8903621556461;14.6221819960590;14.5192877521914;14.5973279759261;14.9476899643253;15.4233504315722;15.9355274033197;16.5116252947223;17.1270957526536;17.7007647391428;18.0859170419278;18.3063513390262;18.3963382467977;18.3039383915083;18.0257824158344;17.5994527413833;17.0890556994325;16.5070300627818;15.9311806007247;15.4197093199240;14.9450348702671;14.5958892051093;14.5192204125840;14.6234912665799;14.8929023380875;15.3101463430400;15.8888685475622;16.5084806025024;10.7623953590755;11.2901846787904;11.7796460325207;12.1321448904323;12.1736536143724;12.0203554539248;11.7401174863064;11.3129760978568;10.7709989288910;10.1715762079976;9.56648675671436;9.01007268222894;8.65191195307112;8.44266556300323;8.29385914269345;8.30718222053088;8.59988135159353;9.03750831470356;9.57725471659727;10.1753661867433;0.603930991971226;1.15454965088465;1.52568456324137;1.76001080986865;1.93327885102392;1.92983454032879;1.61162677627208;1.14288613467936;0.599147559834602;8.68115946043002e-05;-0.598988037848905;-1.14276405297462;-1.61155375579238;-1.92982186320596;-1.93334168349185;-1.76015781483728;-1.52591149560665;-1.15483681806960;-0.604257047252114;-0.000169884159713214;-9.57706352859525;-9.03734942849586;-8.59977372412056;-8.30713768228199;-8.29388420587330;-8.44276020231786;-8.65207245388481;-9.01028509538148;-9.56672835783919;-10.1718277366301;-10.7712422945609;-11.3131875409515;-11.7402768229512;-12.0204500253251;-12.1736778529292;-12.1320937184440;-11.7795292946700;-11.2900189808615;-10.7621990760045;-10.1751629556414;-15.8888846764441;-15.3101529846995;-14.8928992190615;-14.6234752334876;-14.5191922606200;-14.5958545752180;-14.9449893037570;-15.4196495080237;-15.9311059050930;-16.5069451581127;-17.0889691900173;-17.5993721586900;-18.0257142156945;-18.3038872331053;-18.3963088024392;-18.3063460123055;-18.0859310574011;-17.7007895832167;-17.1271236445323;-16.5085043620889;-15.9354043952716;-15.4232446446637;-14.9476104197307;-14.5972837991575;-14.5192821405721;-14.6222128435878;-14.8904221679445;-15.3067216654527;-15.8848610798295;-16.5042587214573;-17.1230057131884;-17.6972384633390;-18.0834106341987;-18.3050781449424;-18.3963827614680;-18.3053183557260;-18.0283544466253;-17.6029503528725;-17.0932402374518;-16.5114970091528;-9.56639354094761;-9.00995499205705;-8.65175684826131;-8.44248655059067;-8.29368121459602;-8.30699446456600;-8.59966472640591;-9.03726249829722;-9.57698694152106;-10.1750852448831;-10.7621062730379;-11.2898861233463;-11.7793442044128;-12.1318702918009;-12.1734601594477;-12.0202151656243;-11.7400072681084;-11.3128838479433;-10.7709160337710;-10.1714933092662;0.599048555629691;1.14281549809409;1.61159139467772;1.92983338482413;1.93331801915188;1.76009758785000;1.52582154290291;1.15473184398399;0.604146013106605;5.84598354854726e-05;-0.604033371937027;-1.15463111243780;-1.52573835803017;-1.76003740875306;-1.93328327155531;-1.92982307406397;-1.61160268038654;-1.14284370408458;-0.599089150620833;-2.25845592714253e-05];
xi=x0;

%x0 = vertcat(cellpoint_out(:,1),cellpoint_out(:,2),cellpoint(:,1),cellpoint(:,2));    %Initial Position of Points
options = odeset('RelTol',1e-5)%,'OutputFcn',@odeplot);                                 %Option Settings for solver

for total_t = first_step:time_step:simulation_length                                  %Total Simulation Time And Step Size
    
    
    interval=interval+1;                                                              %Next Interval
    [t,newpoints] = ode45(@Acini_Function_v1_T1,[start:0.1:total_t],x0,options);             %ODE45 Solver
    
    p = newpoints(size(t,1),:).';                                                     %Initial Positions Of All Points
    new_point(:,interval) = newpoints(size(t,1),:).';                                 %Stores Position Values On Workspace
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Force Calculation on Cell Points From Neighboring Cell Points and Nucleus Points
    
    for s=1:(ncell*npoints)                %Iteration, Points 1 Through (ncell*npoints) of Cell Points
        
        
        if Cell_Point_Number(s)>1 && Cell_Point_Number(s)<20                                                                             %If Cell Point # 1<CP<npoints use this for loop
            d1(s)=sqrt((p(s)-p(s+1))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
            d2(s)=sqrt((p(s)-p(s-1))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
            F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(p(s+1)-p(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(p(s-1)-p(s))/d2(s);                       %Force Calculation Between Cell Points X-Direction
            F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(p(s+(ncell*npoints)+1)-p(s+(ncell*npoints)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
                k_cout(s)*(d2(s)-l0_out)*(p(s+(ncell*npoints)-1)-p(s+(ncell*npoints)))/d2(s);
            
        elseif Cell_Point_Number(s) == 1                                                                                                 %If Cell Point # = 1 use this for loop
            d1(s)=sqrt((p(s)-p(s+1))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
            d2(s)=sqrt((p(s)-p(s+(npoints-1)))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)+(npoints-1)))^2);                             %Distance from Cell Point to Cell Point Right
            F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(p(s+1)-p(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(p(s+(npoints-1))-p(s))/d2(s);             %Force Calculation Between Cell Points X-Direction
            F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(p(s+(ncell*npoints)+1)-p(s+(ncell*npoints)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
                k_cout(s)*(d2(s)-l0_out)*(p(s+(ncell*npoints)+(npoints-1))-p(s+(ncell*npoints)))/d2(s);
            
        elseif Cell_Point_Number(s) == npoints                                                                                           %If Cell Point # = npoints use this for loop
            d1(s)=sqrt((p(s)-p(s-(npoints-1)))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)-(npoints-1)))^2);                             %Distance from Cell Point to Cell Point Left
            d2(s)=sqrt((p(s)-p(s-1))^2+(p(s+(ncell*npoints))-p(s+(ncell*npoints)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
            F_neighx(s)=k_cout(s)*(d1(s)-l0_out)*(p(s-(npoints-1))-p(s))/d1(s)+k_cout(s)*(d2(s)-l0_out)*(p(s-1)-p(s))/d2(s);
            F_neighy(s)=k_cout(s)*(d1(s)-l0_out)*(p(s+(ncell*npoints)-(npoints-1))-p(s+(ncell*npoints)))/d1(s)+...                       %Force Calculation Between Cell Points Y-Direction
                k_cout(s)*(d2(s)-l0_out)*(p(s+(ncell*npoints)-1)-p(s+(ncell*npoints)))/d2(s);
            
        end
        
        d0(s)=sqrt((p(s)-p(s+(2*ncell*npoints)))^2+(p(s+(ncell*npoints))-p(s+(3*ncell*npoints)))^2);       %Distance From Cell Point to Corresponding Nucleus Point
        
        F_outinx(s)=k_in(s)*(d0(s)-l0_inout)*(p(s+(2*ncell*npoints))-p(s))/d0(s);                          %Force Calculation Between Cell Point and Nucleus Point X-Direction
        F_outiny(s)=k_in(s)*(d0(s)-l0_inout)*(p(s+(3*ncell*npoints))-p(s+(ncell*npoints)))/d0(s);          %Force Calculation Between Cell Point and Nucleus POint Y-Direction
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Initialization
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Force Calculation on Cell to Cell Adhesion
        
        F_adhx=0;                                   %Zeros all Cell to Cell Adhesions Forces
        F_adhy=0;                                   %Zeros all Cell to Cell Adhesions Forces
        
        if attached(s)==0                                                       %If Cell Point Is Not Attached Go Through This Loop
            for w=1:ncell*npoints                                               %Iteration, Points 1 through (ncell*npoints) of Cell Structure
                if c(w)~=c(s) && attached(w) == 0 && attached(s)==0             %If Points Not From Same Cell and Are Not Attached Proceed
                    
                    dist(s,w)=(sqrt((p(s)-p(w))^2+(p(s+ncell*npoints)-p(w+ncell*npoints))^2));                        %Distance From Cell Point To All Other Points of Sorrounding Cell Points
                    if dist(s,w)<thresh_max  %&& dist(s,w)>thresh_min                                                   %If Cell Points Within Distance Threshold Proceed
                        attached(s)=w;                                                                                %Cell Point S Becomes Attached
                        attached(w)=s;                                                                                %Cell Point W Becomes Attached
                        F_adhx=F_adhx+k_cell*(dist(s,w)-l_bond)*(p(w)-p(s))/dist(s,w);                                %Adhesion Forces in X-direction
                        F_adhy=F_adhy+k_cell*(dist(s,w)-l_bond)*(p(w+ncell*npoints)-p(s+ncell*npoints))/dist(s,w);    %Adhesion Forces in Y-direction
                        
                    else
                        F_adhx=F_adhx+0;                                        %Adhesion Force Unchanged
                        F_adhy=F_adhy+0;                                        %Adhesion Force Unchanged
                    end
                    
                else
                    F_adhx=F_adhx+0;                                        %Adhesion Force Unchanged
                    F_adhy=F_adhy+0;
                end
            end
            %if c(w)~=c(s) && attached(w) == s && attached(s)==w
            %dist(s,w)=(sqrt((i(s)-i(w))^2+(i(s+ncell*npoints)-i(w+ncell*npoints))^2));
            %if dist(s,w)<new_cell_thresh && dist(s,w)>new_cell_thresh
            %add new cell
            %alter the npoints ncells in correct order
            %end
            
            
        elseif attached(s) ~= 0
            w = attached(s);
            if F_adh_out(s) > 4
                k_in(s) = (F_adh_out(s)/(est_k + F_adh_out(s))) * k_in_new(s) ;
            end
            dist(s,w)=(sqrt((p(s)-p(w))^2+(p(s+ncell*npoints)-p(w+ncell*npoints))^2));                        %Distance From Cell Point To All Other Points of Sorrounding Cell Points
            if dist(s,w)<thresh_max  %&& dist(s,w)>thresh_min                                                 %If Cell Points Within Distance Threshold Proceed
                
                F_adhx=F_adhx+k_cell*(dist(s,w)-l_bond)*(p(w)-p(s))/dist(s,w);                                %Adhesion Forces in X-direction
                F_adhy=F_adhy+k_cell*(dist(s,w)-l_bond)*(p(w+ncell*npoints)-p(s+ncell*npoints))/dist(s,w);    %Adhesion Forces in Y-direction
                
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
        F_adh_out(s) = sqrt((F_adhx_out(s))^2+(F_adhy_out(s))^2);
        if attached(s)~=0                                                   %If Attached Proceed
            F_adhx_out(attached(s)) = -F_adhx;                              %Opposite Reaction Force X-Direction
            F_adhy_out(attached(s)) = -F_adhy;                              %Opposite Reaction Force X-Direction
        end
        Fx_adh_matrix (s, interval) = F_adhx_out (s);                      %Stores Force Values On Workspace
        Fy_adh_matrix (s, interval) = F_adhy_out (s);                      %Stores Force Values On Workspace           
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         if attached(s)~=0                                                   %If Attached Proceed
        %             F_adhx_out(attached(s)) = -F_adhx;                              %Opposite Reaction Force X-Direction
        %             F_adhy_out(attached(s)) = -F_adhy;                              %Opposite Reaction Force X-Direction
        %         end
        %
        if attached(s) ~= 0
            
            F_extx (s) = 0;                                                 %Zero Force In X-Direction
            F_exty (s) = 0;                                                 %Zero Force In Y-Direction
            
        elseif attached(s) == 0
            
            r_x = mean(p(((2*ncell*npoints)+1):(3*ncell*npoints)));                              %Mean X-Position of Acini Structure
            r_y = mean(p(((3*ncell*npoints)+1):(4*ncell*npoints)));                              %Mean Y-Position of Acini Structure
            dist_center(s) = sqrt((p(s)-r_x)^2+(p(s+ncell*npoints)-r_y)^2);                      %Distance of Cell Point From Center of Acini Structure
            
            if c(s) == 1
                
                nr_x = mean(p(((2*ncell*npoints)+1):(2*ncell*npoints+npoints)));
                nr_y = mean(p(((3*ncell*npoints)+1):(3*ncell*npoints+npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 2
                
                nr_x = mean(p(((2*ncell*npoints)+npoints+1):(2*ncell*npoints+2*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+npoints+1):(3*ncell*npoints+2*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 3
                
                nr_x = mean(p(((2*ncell*npoints)+2*npoints+1):(2*ncell*npoints+3*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+2*npoints+1):(3*ncell*npoints+3*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 4
                
                nr_x = mean(p(((2*ncell*npoints)+3*npoints+1):(2*ncell*npoints+4*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+3*npoints+1):(3*ncell*npoints+4*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 5
                
                nr_x = mean(p(((2*ncell*npoints)+4*npoints+1):(2*ncell*npoints+5*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+4*npoints+1):(3*ncell*npoints+5*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 6
                
                nr_x = mean(p(((2*ncell*npoints)+5*npoints+1):(2*ncell*npoints+6*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+5*npoints+1):(3*ncell*npoints+6*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 7
                
                nr_x = mean(p(((2*ncell*npoints)+6*npoints+1):(2*ncell*npoints+7*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+6*npoints+1):(3*ncell*npoints+7*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 8
                
                nr_x = mean(p(((2*ncell*npoints)+7*npoints+1):(2*ncell*npoints+8*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+7*npoints+1):(3*ncell*npoints+8*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 9
                
                nr_x = mean(p(((2*ncell*npoints)+8*npoints+1):(2*ncell*npoints+9*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+8*npoints+1):(3*ncell*npoints+9*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            elseif c(s) == 10
                
                nr_x = mean(p(((2*ncell*npoints)+9*npoints+1):(2*ncell*npoints+10*npoints)));
                nr_y = mean(p(((3*ncell*npoints)+9*npoints+1):(3*ncell*npoints+10*npoints)));
                dist_n(s) = sqrt((p(s)-nr_x)^2+(p(s+ncell*npoints)-nr_y)^2);
                
            end
            
            %dist_nr = 15;  %sqrt((nr_x-r_x)^2+(nr_y-r_y)^2);
            
            if F_external(s) ~= 0  && dist_center(s) > dist_nr                                    %External Force Generated
                F_extx (s) = F_external(s)*(p(s)-nr_x)/dist_n(s);
                F_exty (s) = F_external(s)*(p(s+(ncell*npoints))-nr_y)/dist_n(s);
                
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
        
        Fx_out_matrix (s, interval) = Fx_out (s);                      %Stores Force Values On Workspace
        Fy_out_matrix (s, interval) = Fy_out (s);                      %Stores Force Values On Workspace
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Force Calculation On Nucleus Points From Neighboring Nucleus Points and Cell Points
    
    for q=1:(ncell*npoints) %Iteration, Points 1 Through (ncell*npoints) of Nucleus
        
        if  Cell_Point_Number(q) > 1 && Cell_Point_Number(q) < 20                                                               %If Nucleus Point # 1 < CP < npoints Use This For Loop
            
            d1_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)+1))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)+1))^2);       %Distance from Cell Point To Nucleus Point Left
            d2_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)-1))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)-1))^2);       %Distance from Cell Point To Nucleus Point Right
            F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)+1)-p(q+(2*ncell*npoints)))/d1_in(q)+...                                      %Force Calculation Between Nucleus Points X-Direction
                k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)-1)-p(q+(2*ncell*npoints)))/d2_in(q);
            F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)+1)-p(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)+...      %Force Calculation Between Nucleus Points Y-Direction
                k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)-1)-p(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
            
        elseif Cell_Point_Number(q) == 1                                                                                        %If Nucleus Point # = 1 Use This For Loop
            
            d1_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)+1))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)+1))^2);                     %Distance from Cell Point To Nucleus Point Left
            d2_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)+(npoints-1)))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)+(npoints-1)))^2); %Distance from Cell Point To Nucleus Point Right
            F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)+1)-p(q+(2*ncell*npoints)))/d1_in(q)+...                                                    %Force Calculation Between Nucleus Points X-Direction
                k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(npoints-1))-p(q+(2*ncell*npoints)))/d2_in(q);
            F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)+1)-p(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)+...                    %Force Calculation Between Nucleus Points Y-Direction
                k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)+(npoints-1))-p(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
            
        elseif Cell_Point_Number(q) == npoints                                                                                  %If Nucleus Point # = npoints Use This For Loop
            
            d1_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)-(npoints-1)))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)-(npoints-1)))^2); %Distance from Cell Point To Nucleus Point Left
            d2_in(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)-1))^2+(p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)-1))^2);                     %Distance from Cell Point To Nucleus Point Right
            F_inneighx(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)-(npoints-1))-p(q+(2*ncell*npoints)))/d1_in(q)+...                                          %Force Calculation Between Nucleus Points X-Direction
                k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)-1)-p(q+(2*ncell*npoints)))/d2_in(q);
            F_inneighy(q)=k_cin(q+(2*ncell*npoints))*(d1_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)-(npoints-1))-p(q+(2*ncell*npoints)+(ncell*npoints)))/d1_in(q)...           %Force Calculation Between Nucleus Points Y-Direction
                +k_cin(q+(2*ncell*npoints))*(d2_in(q)-l0_in)*(p(q+(2*ncell*npoints)+(ncell*npoints)-1)-p(q+(2*ncell*npoints)+(ncell*npoints)))/d2_in(q);
            
        end
        
        d0(q)=sqrt((p(q)-p(q+(2*ncell*npoints)))^2+(p(q+(ncell*npoints))-p(q+(3*ncell*npoints)))^2);           %Distance From Cell Point to Corresponding Nucleus Point
        
        F_inoutx(q)=k_in(q)*(d0(q)-l0_inout)*(p(q)-p(q+(2*ncell*npoints)))/d0(q);                              %Force Calculation Between Cell Point and Nucleus Point X-Direction
        F_inouty(q)=k_in(q)*(d0(q)-l0_inout)*(p(q+(ncell*npoints))-p(q+(3*ncell*npoints)))/d0(q);              %Force Calculation Between Cell Point and Nucleus Point Y-Direction
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Force Calculation On Internal Nucleus
        
        
        if  Cell_Point_Number(q) >= 1 && Cell_Point_Number(q) <= 10                                                         %If Nucleus Point Is >=1 & <=10 Use This For Loop
            
            d_nu(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)+10))^2+...                                           %Distance calculation Between Nucleus Points
                (p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)+10))^2);
            
            F_nx(q)=k_nu(q+(2*ncell*npoints))*(d_nu(q)-l0_nu)*(p(q+(2*ncell*npoints)+10)-p(q+(2*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
            F_ny(q)=k_nu(q+(3*ncell*npoints))*(d_nu(q)-l0_nu)*(p(q+(3*ncell*npoints)+10)-p(q+(3*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
            
        elseif Cell_Point_Number(q) >= 11 && Cell_Point_Number(q) <= 20                                                     %If Nucleus Point Is >=11 & <=20 Use This For Loop
            
            d_nu(q)=sqrt((p(q+(2*ncell*npoints))-p(q+(2*ncell*npoints)-10))^2+...                                           %Distance calculation Between Nucleus Points
                (p(q+(2*ncell*npoints)+(ncell*npoints))-p(q+(2*ncell*npoints)+(ncell*npoints)-10))^2);
            
            F_nx(q)=k_nu(q+(2*ncell*npoints))*(d_nu(q)-l0_nu)*(p(q+(2*ncell*npoints)-10)-p(q+(2*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
            F_ny(q)=k_nu(q+(3*ncell*npoints))*(d_nu(q)-l0_nu)*(p(q+(3*ncell*npoints)-10)-p(q+(3*ncell*npoints)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
            
        end
        
        F_inadhx=0;       %Zero Force In X-Direction
        F_inadhy=0;       %Zero Force In Y-Direction
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Total Nucleus Force Calculation
        
        Fx_in(q)= F_inneighx(q)+F_inoutx(q)+F_inadhx+F_intx(q)+F_nx(q);    %Sum Of Forces In X-Direction
        Fy_in(q)= F_inneighy(q)+F_inouty(q)+F_inadhy+F_inty(q)+F_ny(q);    %Sum Of Forces In X-Direction
        
        Fx_in_matrix(q,interval) = Fx_in(q);                               %Stores Force Values On Workspace
        Fy_in_matrix(q,interval) = Fy_in(q);                               %Stores Force Values On Workspace
        
    end
    
    F_in_mag(:,interval)=sqrt(Fx_in_matrix(:,interval).^2+Fy_in_matrix(:,interval).^2);     %Stores Force Values Magnitude On Workspace
    F_out_mag(:,interval)=sqrt(Fx_out_matrix(:,interval).^2+Fy_out_matrix(:,interval).^2);  %Stores Force Values Magnitude On Workspace
    F_adh_mag(:,interval)=sqrt(Fx_adh_matrix(:,interval).^2+Fy_adh_matrix(:,interval).^2);  %Stores Force Values Magnitude On Workspace    
    forces(:,interval) = F_external(:,1).';                                                 %Stores Force Values Magnitude On Workspace
    
    start = total_t;                                                                        %Sets New Simulation Start Time
    x0 = p.' ;                                                                              %Sets New Simulation Start Positions
    force_number = force_number + 1;                                                        %Sets New Simulation Random Force Vector
    
    point_set_1 = [];
    point_set_2 = [];
    point_set_3 = [];
    point_set_4 = [];
    point_set_5 = [];
    point_set_6 = [];
    point_set_7 = [];
    point_set_8 = [];
    point_set_9 = [];
    point_set_10 = [];
    
    point_set = [];
    new_forcepoint = zeros(1,200);
    for sp=1:(ncell*npoints)
        if attached(sp) == 0
            
            r_x = mean(p(((2*ncell*npoints)+1):(3*ncell*npoints)));                              %Mean X-Position of Acini Structure
            r_y = mean(p(((3*ncell*npoints)+1):(4*ncell*npoints)));                              %Mean Y-Position of Acini Structure
            dist_center(sp) = sqrt((p(sp)-r_x)^2+(p(sp+ncell*npoints)-r_y)^2);
            
            if dist_center(sp) > dist_nr
                
                new_forcepoint(sp)=sp;
                
            else
                
                new_forcepoint(sp)=0;
            end
        else
            new_forcepoint(sp)=0;
        end
    end
    
    for sp=1:(ncell*npoints)
        
        if c(sp)==1 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_1 = [point_set_1 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==2 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_2 = [point_set_2 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==3 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_3 = [point_set_3 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==4 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_4 = [point_set_4 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==5 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_5 = [point_set_5 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==6 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_6 = [point_set_6 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==7 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_7 = [point_set_7 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==8 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_8 = [point_set_8 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==9 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_9 = [point_set_9 new_forcepoint(sp)]  ;
            
        elseif   c(sp)==10 && new_forcepoint(sp)~=0 && attached(sp) == 0
            
            point_set_10 = [point_set_10 new_forcepoint(sp)] ;
            
            
        end
    end
    
    
%     for sig = 1:(ncell*npoints)
%         if c(sig)==1
%             if F_in_mag(sig,interval) > 1.5
%                 
%                 v(1:20) = 10;
%            
%                 ECM1 = 15;
%                 
%                 k_cin(1:20) = 14;
%                 k_cout(1:20) = 5;
%                 k_nu(1:20) = 14;
%                 
%             end
%             
%         elseif c(sig)==2
%             if F_in_mag(sig,interval) > 1.5
%                 
%                 v(21:40) = 10;
%                 
%                 ECM2 = 15;
%                 
%                 k_cin(21:40) = 14;
%                 k_cout(21:40) = 5;
%                 k_nu(21:40) = 14;
%                 
%             end
%             
%         elseif c(sig)==3
%              if F_in_mag(sig,interval) > 1.5
%                 
%                 v(41:60) = 10;
%                 
%                 ECM3 = 15;
%                 
%                 k_cin(41:60) = 14;
%                 k_cout(41:60) = 5;
%                 k_nu(41:60) = 14;
%                 
%                             
%             end
%             
%         elseif c(sig)==4
%              if F_in_mag(sig,interval) > 1.5
%                 
%                 v(61:80) = 10;
%                 
%                 ECM4 = 15;
%                 
%                 k_cin(61:80) = 14;
%                 k_cout(61:80) = 5;
%                 k_nu(61:80) = 14;
%                 
%             
%             end
%             
%         elseif c(sig)==5
%            if F_in_mag(sig,interval) > 1.5
%                 
%                 v(81:100) = 10;
%                 
%                 ECM5 = 15;
%                 
%                 k_cin(81:100) = 14;
%                 k_cout(81:100) = 5;
%                 k_nu(81:100) = 14;
%                 
%             
%             end
%             
%             
%         elseif c(sig)==6
%               if F_in_mag(sig,interval) > 1.5
%                 
%                 v(101:120) = 10;
%                 
%                 ECM6 = 15;
%                 
%                 k_cin(101:120) = 14;
%                 k_cout(101:120) = 5;
%                 k_nu(101:120) = 14;
%                 
%             
%             end
%             
%         elseif c(sig)==7
%              if F_in_mag(sig,interval) > 1.5
%                 
%                 v(121:140) = 10;
%                 
%                 ECM7 = 15;
%                 
%                 k_cin(121:140) = 14;
%                 k_cout(121:140) = 5;
%                 k_nu(121:140) = 14;
%                 
%            
%             end
%             
%         elseif c(sig)==8
%               if F_in_mag(sig,interval) > 1.5
%                 
%                 v(141:160) = 10;
%                 
%                 ECM8 = 15;
%                 
%                 k_cin(141:160) = 14;
%                 k_cout(141:160) = 5;
%                 k_nu(141:160) = 14;
%                 
%             
%             end
%             
%         elseif c(sig)==9
%                if F_in_mag(sig,interval) > 1.5
%                 
%                 v(161:180) = 10;
%                 
%                 ECM9 = 15;
%                 
%                 k_cin(161:180) = 14;
%                 k_cout(161:180) = 5;
%                 k_nu(161:180) = 14;
%                 
%                             
%                end
%             
%         elseif c(sig)==10
%              if F_in_mag(sig,interval) > 1.5
%                 
%                 v(181:200) = 10;
%                 
%                 ECM10 = 15;
%                 
%                 k_cin(181:200) = 14;
%                 k_cout(181:200) = 5;
%                 k_nu(181:200) = 14;
%                 
%            
%             end
%         end
%     end
    
    
    ECM = {ECM1,ECM2,ECM3,ECM4,ECM5,ECM6,ECM7,ECM8,ECM9,ECM10};
    point_set = {point_set_1,point_set_2,point_set_3,point_set_4,point_set_5,point_set_6,point_set_7,point_set_8,point_set_9,point_set_10};
    
    F_external  = Acini_ForceGen_v2_T1(ECM,point_set,SD);
    interval
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots

figure
scatter(xi(1:ncell*npoints),xi(ncell*npoints+1:2*ncell*npoints),[],zeros(1,200),'*')
axis equal
title('Original')
colorbar
ylabel(colorbar, 'Stress nN')
hold on
scatter(xi(2*ncell*npoints+1:3*ncell*npoints),xi(3*ncell*npoints+1:4*ncell*npoints),[],zeros(1,200),'*')
hold off

figure
scatter(newpoints(size(t,1),1:ncell*npoints),newpoints(size(t,1),(ncell*npoints+1):2*ncell*npoints),[],F_out_mag(:,end),'*')
axis equal
title('Final')
colorbar
ylabel(colorbar, 'Stress nN')
hold on
scatter(newpoints(size(t,1),(2*ncell*npoints+1):3*ncell*npoints),newpoints(size(t,1),(3*ncell*npoints+1):4*ncell*npoints),[],F_in_mag(:,end),'*')
hold off

figure
scatter(new_point(1:ncell*npoints,(simulation_length/time_step)/2),new_point((ncell*npoints+1):2*ncell*npoints,(simulation_length/time_step)/2),[],F_out_mag(:,(simulation_length/time_step)/2),'*')
axis equal
title('Original & Midpoint')
colorbar
ylabel(colorbar, 'Stress nN')
hold on
scatter(new_point((2*ncell*npoints+1):3*ncell*npoints,(simulation_length/time_step)/2),new_point((3*ncell*npoints+1):4*ncell*npoints,(simulation_length/time_step)/2),[],F_in_mag(:,(simulation_length/time_step)/2),'*')
hold on
scatter(xi(1:ncell*npoints),xi(ncell*npoints+1:2*ncell*npoints),[],zeros(1,200),'o')
hold on
scatter(xi(2*ncell*npoints+1:3*ncell*npoints),xi(3*ncell*npoints+1:4*ncell*npoints),[],zeros(1,200),'o')
hold off

figure
scatter(newpoints(size(t,1),1:ncell*npoints),newpoints(size(t,1),(ncell*npoints+1):2*ncell*npoints),[],F_out_mag(:,end),'*')
axis equal
title('Original & Final')
colorbar
ylabel(colorbar, 'Stress nN')
hold on
scatter(newpoints(size(t,1),(2*ncell*npoints+1):3*ncell*npoints),newpoints(size(t,1),(3*ncell*npoints+1):4*ncell*npoints),[],F_in_mag(:,end),'*')
hold on
scatter(xi(1:ncell*npoints),xi(ncell*npoints+1:2*ncell*npoints),[],zeros(1,200),'o')
hold on
scatter(xi(2*ncell*npoints+1:3*ncell*npoints),xi(3*ncell*npoints+1:4*ncell*npoints),[],zeros(1,200),'o')
hold off

