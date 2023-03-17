global  Param Qsteadystate qfsteadystate X0 Q hc row ros Vsd0 ph dpdt2 dardt2 dpdt1 dardt1 phin Vsd02 Param2 Qsteadystate2 qfsteadystate2 X02 hf Param3 Ts_out1_init_PB Ts_out2_init_PB Attemp1_PB Ts_out3_init_PB Ts_out4_init_PB Attemp2_PB  Ts_out5_init_PB  Ts_out6_init_PB Param4 Ts_out1_init_RB Ts_out2_init_RB Attemp1_RB Ts_out3_init_RB Ts_out4_init_RB Attemp2_RB  Ts_out5_init_RB  Ts_out6_init_RB Attempflow_init Attempflow_SS Attempflow_SSRB Param5 TFG_out6_initPB TFG_out5_initPB TFG_out4_initPB TFG_out3_initPB TFG_out2_initPB TFG_out1_initPB FGFlow_SS FG_flow_RB Param6 TFG_out6_initRB TFG_out5_initRB TFG_out4_initRB TFG_out3_initRB TFG_out2_initRB TFG_out1_initRB


warning off



phin=6.21000008802327;

%%%%%%%%%%%%%%%Param for RB Boiler%%%%%%%

hf=9e5;
dpdt1  = 0;
dardt1 = 0;

Vsd0 = 2.105229086328673e+000;

qfsteadystate =73.2301;

Qsteadystate  = 1.374070083420811e+008;
X0 =[106.953303000225
    6700000.56080471
    0.0429635897773415
     1.61942807169209];
% Adc (m2)
Param(1) = 0.33;
% Vd (m3)
Param(2) = 20.4;
% Vdc (m3)
Param(3) = 26.9;
% Ad (m2)
Param(4) = 7.64;
% beta (dimensionless)
Param(5) = 0.3;
% k (dimensionless)       
Param(6) = 25;
% g (m/s2)      
Param(7) = 9.81;
% Td (s)     
Param(8) = 1;
% cp (J/K.kg)
Param(9) = 550;
% mt (kg)      
Param(10) = 480185.6;
% mr (kg)     
Param(11) = 345156.9;
% Vr (m3)
Param(12) = 108.5;
% md (kg)
Param(13) = 50909.1;


%%%%%%%%%%%%%%%Param2 for PB Boiler%%%%%%%
hf2=9e5;

dpdt2  = 0;
dardt2 = 0;

Vsd02 = 1.622713827338826e+000;

qfsteadystate2 =31.5677;
Qsteadystate2  =5.9194e+007;
X02 =[39.2837325596915
    6800297.15097339
    0.0336035188446647
     1.24093404723588];

% Adc (m2)
Param2(1) = 0.33;
% Vd (m3)
Param2(2) = 14.8;

% Vdc (m3)
Param2(3) = 6.75;
% Ad (m2)
Param2(4) = 5.6;
% beta (dimensionless)
Param2(5) = 0.3;
% k (dimensionless)       
Param2(6) = 25;
% g (m/s2)      
Param2(7) = 9.81;
% Td (s)     
Param2(8) = 10/5.5;
% cp (J/K.kg)
Param2(9) = 550;
% mt (kg)      
Param2(10) = 198000;
% mr (kg)     
Param2(11) = 140000;
% Vr (m3)
Param2(12) = 40;
% md (kg)
Param2(13) = 38000;



%%%%%%%%%%%%%%%Param3 for PB Superheater%%%%%%%
%Tube diameter, r_tube (m) Taler
Param3 (1)= 2.1e-2;

%Tube length, L_tube (m) per compartment Taler
Param3 (2)= 0.743;

%T_FG_out1
Param3 (3) = 471; %oC
%T_FG_out2
Param3 (4) = 534; %oC
%Mass of superheater tube per compartment (kg)
Param3 (5)= 4.55E3; %Taler

%attemp1 radius, r_attemp1 (m)
Param3 (6)= 0.17; %Cho Baekhyun

%attemp2 length, L_tube (m)
Param3 (7)= 2.6; %Cho Baekhyun


%Water temp attemporator (degC)
Param3 (8)= 100; %Zima.W

%Attemporator Mass (kg) Cho Baekhyun
Param3 (9)= 245.85;

%Attemporator Pressure(bar)
Param3 (10) =115; %Zima W

%number of tube SH 
Param3(11)=148; %Taler

%T_FG_out3
Param3 (12) = 598; %oC

%T_FG_out4
Param3 (13) = 661; %oC

%T_FG_out5
Param3 (14) = 724; %oC

%T_FG_out6
Param3 (15) = 790; %oC

%Cp_Metal
Param3 (16)= 0.466; %kJ/kgoC https://gchem.cm.utexas.edu/data/section2.php?target=heat-capacities.php

%Attemporator flow initial (kg/hr)
Attempflow_init =0;
%Atemperator Flow Steady State (kg/hr)
Attempflow_SS = 11355;

%Initial Condition of Ts_out1
Ts_out1_init_PB = 303.33953 %degC

Ts_out2_init_PB = 334.14047;%degC

Attemp1_PB = 334.14047;

Ts_out3_init_PB = 375.4672; %degC

Ts_out4_init_PB = 426.4464; %degC

Attemp2_PB =426.4464;

Ts_out5_init_PB = 485.4202%degC

Ts_out6_init_PB = 550.8115 %degC


%%%%%%%%%%%%%%%Param3 for RB Superheater%%%%%%%

%Tube diameter, r_tube (m)
Param4 (1)= 2.1e-2;

%Tube length, L_tube (m) per compartment
Param4 (2)= 0.743;

%T_FG_out1
Param4 (3) = 471; %oC
%T_FG_out2
Param4 (4) = 534; %oC
%Mass of superheater tube per compartment (kg)
Param4 (5)= 1.06E4; %Taler

%attemp1 radius, r_attemp1 (m)
Param4 (6)= 0.17; %Taler

%attemp2 length, L_tube (m)
Param4 (7)= 2.6; %Taler


%Water temp attemporator (degC)
Param4 (8)= 150; %Zima.W

%Attemporator Mass (kg) Cho, Baekhyun
Param4 (9)= 245.85;

%Attemporator Pressure(bar) Zima W
Param4 (10) =115; %Zima W

%number of tube SH 
Param4(11)=330; %Taler

%T_FG_out3
Param4 (12) = 598; %oC

%T_FG_out4
Param4 (13) = 661; %oC

%T_FG_out5
Param4 (14) = 724; %oC

%T_FG_out6
Param4 (15) = 790; %oC

%Cp_Metal
Param4(16)= 0.466; %kJ/kgoC https://gchem.cm.utexas.edu/data/section2.php?target=heat-capacities.php
%Atemperator Flow Steady State (kg/hr)
Attempflow_SSRB = 2.3*11355;

%Initial Condition of Ts_out1
Ts_out1_init_RB = 302.2065815; %degC

Ts_out2_init_RB =332.70172955;%degC

Attemp1_RB= 332.70172955;

Ts_out3_init_RB = 373.4655693; %degC

Ts_out4_init_RB =423.6721137;  ;%degC


Attemp2_RB= 423.6721137;
Ts_out5_init_RB = 481.745302 ;%degC

Ts_out6_init_RB = 546.1691598 %degC


%%%%%%%%%%%%%%%Param5 for PB Superheater FG Dynamic%%%%%%%
%Tube diameter, r_tube (m) Taler
Param5 (1)= 2.1e-2;

%Tube length, L_tube (m) per compartment Taler
Param5 (2)= 0.743;

%Mass of superheater tube per compartment (kg)
Param5 (3)= 4.55E3;  %Taler

%number of tube SH 
Param5(4)=148; %Taler

%Cp of FG assummed it's constant
Param5 (5)= 1.144; %kJ/kgK Ordys


%Cp_Metal
Param5 (6)= 0.466; %kJ/kgoC https://gchem.cm.utexas.edu/data/section2.php?target=heat-capacities.php

%Area of convection flue gas excluded Tube per Tube(Taler) , (m2)
Param5 (7) = 1.73E-2;

%Density of flue gas (assumed it's constant) kg/m3
Param5 (8) = 0.4 ;


TFG_out6_initPB= 764.5389; %oC
TFG_out5_initPB= 684.8106; %oC
TFG_out4_initPB= 611.0105; %oC
TFG_out3_initPB= 543.7293; %oC
TFG_out2_initPB= 483.8618; %oC
TFG_out1_initPB= 432.2968; %oC

FGFlow_SS = 182700; %back calculation of Steady State
FG_flow_RB =2.3*FGFlow_SS;

%%%%%%%%%%%%%%%Param5 for RB Superheater FG Dynamic%%%%%%%
%Tube diameter, r_tube (m) Taler
Param6 (1)= 2.1e-2;

%Tube length, L_tube (m) per compartment Taler
Param6 (2)= 0.743;

%Mass of superheater tube per compartment (kg)
Param6 (3)= 1.06E4; %Taler

%number of tube SH 
Param6(4)=330; %Taler

%Cp of FG assummed it's constant
Param6 (5)= 1.144; %kJ/kgK Ordys

%Cp_Metal
Param6(6)= 0.466; %kJ/kgoC https://gchem.cm.utexas.edu/data/section2.php?target=heat-capacities.php

%Area of convection flue gas excluded Tube per Tube(Taler) , (m2)
Param6 (7) = 1.73E-2;


%Density of flue gas (assumed it's constant) kg/m3
Param6 (8) = 0.4 ;


TFG_out6_initRB= 765.12376925; %oC
TFG_out5_initRB= 685.9654778; %oC
TFG_out4_initRB= 612.6928895; %oC
TFG_out3_initRB= 545.8638828; %oC
TFG_out2_initRB= 486.31627357; %oC
TFG_out1_initRB= 434.8845734; %oC


