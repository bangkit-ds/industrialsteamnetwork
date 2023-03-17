function [sys,x0,str,tss]=PBsuperheater2(t,x,u,flag,Param3,X_ss)

global  Param3 
switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
        % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param3);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param3);

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

global  Param3 

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 1;     % outputs
sizes.NumInputs      = 5;     % inputs 
sizes.DirFeedthrough = 1;     % System is non-causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the continuous states.
str = [];                     % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].



% ******************************************************************************************************************************
%  Outputs
% ******************************************************************************************************************************


function sys = mdlOutputs(t,x,u,Param3);

global  Param3 
% External Parameters:
r_tube = Param3(1);
L_tube = Param3(2); %tube length
%T_FG_out2 =Param3 (4);
M_tube = Param3 (5);
n_tube=Param3 (11);
Cp_M = Param3(16);


% States
Ts_out2= x(1);

% Inputs
PB_p = u(1)/14.7; %bar
steam_flow= u(2)/3600; %kg/s
U = u(3); %220 W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_in2= u(4);
T_FG_out2 = u(5);

%properties 
rho_SHPB_in= XSteam ('rho_pT',PB_p,Ts_in2); %density steam (in) from previous compartment
rho_SHPB_out= XSteam ('rho_pT',PB_p, Ts_out2); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',PB_p,Ts_in2) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',PB_p,Ts_out2) ; %Cp steam T_out


%Update Derivative


% Outputs:
Cp_s= (Cp_s_in+Cp_s_out)/2;%Cp steam average
rho_s= (rho_SHPB_in+rho_SHPB_out)/2; %density steam average


% sys(1)is dT_out/dt 

sys(1) =  Ts_out2;

%Add other outputs


% ******************************************************************************************************************************
% Derivatives
% ******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param3)

global  Param3 

% External Parameters:
r_tube = Param3(1);
L_tube = Param3(2); %tube length
%T_FG_out2 =Param3 (4);
M_tube = Param3 (5);
n_tube=Param3 (11);
Cp_M = Param3(16);


% States
Ts_out2= x(1)

% Inputs
PB_p = u(1)/14.7; %bar
steam_flow= u(2)/3600; %kg/s
U = u(3); %220 W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_in2= u(4);
T_FG_out2 = u(5);

%properties 
rho_SHPB_in= XSteam ('rho_pT',PB_p,Ts_in2); %density steam (in) from previous compartment
rho_SHPB_out= XSteam ('rho_pT',PB_p, Ts_out2); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',PB_p,Ts_in2) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',PB_p,Ts_out2) ; %Cp steam T_out

% Outputs:
Cp_s= (Cp_s_in+Cp_s_out)/2%Cp steam average
rho_s= (rho_SHPB_in+rho_SHPB_out)/2 %density steam average

% Derivatives:
% sys(1)is dT_out/dt 

sys(1) =(steam_flow *(Cp_s*Ts_in2-Cp_s*Ts_out2) + U*(2*pi*r_tube*L_tube*n_tube)*(T_FG_out2-Ts_out2))/(rho_s*(pi*r_tube^2*L_tube*n_tube)*Cp_s+M_tube*Cp_M);





