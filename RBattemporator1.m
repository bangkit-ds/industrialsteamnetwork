function [sys,x0,str,tss]=RBattemporator1(t,x,u,flag,Param4,X_ss)

global  Param4 
switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
        % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param4);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param4);

otherwise,
   sys = [];
end

% ******************************************
% Sub-routines or Functions
% ******************************************

% ******************************************
% Initialization
% ******************************************
function [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);

global  Param4

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 2;     % outputs
sizes.NumInputs      = 4;     % inputs 
sizes.DirFeedthrough = 1;     % System is non-causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the continuous states.
str = [];                     % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].



% ******************************************************************************************************************************
%  Outputs
% ******************************************************************************************************************************


function sys = mdlOutputs(t,x,u,Param4);

global  Param4
% External Parameters:
r_attempt1 = Param4(6);
L_attemp1 = Param4(7); %tube length
T_water= Param4(8);
M_attemporator = Param4 (9);
P_water = Param4 (10);
Cp_M = Param4 (16);


% States
Ts_out_attemp1= x(1);
% steam_flow_out_attemp1 = x(2);

% Inputs
RB_p = u(1)/14.7; %bar
steam_flow= u(2)/3600; %kg/s
Ts_in_attemp1= u(3); %oC
% water_attemp1_flow= max([u(4)/3600 0]);
water_attemp1_flow = u(4)/3600; %kg/s 

%properties 
rho_SHRB_in= XSteam ('rho_pT',RB_p,Ts_in_attemp1); %density steam (in) from previous compartment
rho_SHRB_out= XSteam ('rho_pT',RB_p, Ts_out_attemp1); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',RB_p,Ts_in_attemp1) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',RB_p,Ts_out_attemp1) ; %Cp steam T_out
Cp_attemp_in=XSteam('Cp_pT',P_water,T_water); %Cp attemp

% Auxiliary Variables
steam_flow_out_attemp1 = steam_flow + water_attemp1_flow; %kg/s

% Outputs:
% states
sys(1) =  Ts_out_attemp1;
sys(2)= steam_flow_out_attemp1;

%auxilary

% ******************************************************************************************************************************
% Derivatives
% ******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param4)

global  Param4 


% External Parameters:
r_attempt1 = Param4(6);
L_attemp1 = Param4(7); %tube length
T_water= Param4(8);
M_attemporator = Param4 (9);
P_water = Param4 (10);
Cp_M = Param4 (16);


% States
Ts_out_attemp1= x(1);
% steam_flow_out_attemp1 = x(2);

% Inputs
RB_p = u(1)/14.7; %bar
steam_flow= u(2)/3600; %kg/s
Ts_in_attemp1= u(3); %oC
% water_attemp1_flow= max([u(4)/3600 0]);
water_attemp1_flow = u(4)/3600; %kg/s 

%properties 
rho_SHRB_in= XSteam ('rho_pT',RB_p,Ts_in_attemp1); %density steam (in) from previous compartment
rho_SHRB_out= XSteam ('rho_pT',RB_p, Ts_out_attemp1); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',RB_p,Ts_in_attemp1) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',RB_p,Ts_out_attemp1) ; %Cp steam T_out
Cp_attemp_in=XSteam('Cp_pT',P_water,T_water); %Cp attemp

Cp_s = (Cp_s_in+Cp_s_out)/2;
rho_s= (rho_SHRB_in+rho_SHRB_out)/2; %density steam average

% Auxiliary Variables
steam_flow_out_attemp1 = steam_flow + water_attemp1_flow; %kg/s

% Derivatives:
sys (1)= (steam_flow*Cp_s*Ts_in_attemp1+ water_attemp1_flow*Cp_attemp_in*T_water-steam_flow_out_attemp1*Cp_s*Ts_out_attemp1)/(rho_s*(pi*r_attempt1^2*L_attemp1)*Cp_s+M_attemporator*Cp_M);
% sys (1)= (steam_flow*(Cp_s*Ts_in_attemp1-Cp_s*Ts_out_attemp1) + water_attemp1_flow*(Cp_attemp_in*T_water-Cp_s*Ts_out_attemp1))/(rho_s*(pi*r_attempt1^2*L_attemp1)*Cp_s+M_attemporator*Cp_M);

% sys (2)= steam_flow + water_attemp1_flow - steam_flow_out_attemp1; %kg/s
