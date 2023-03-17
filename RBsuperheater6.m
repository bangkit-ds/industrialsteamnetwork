function [sys,x0,str,tss]=RBsuperheater6(t,x,u,flag,Param4,X_ss)

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


function sys = mdlOutputs(t,x,u,Param4);

global  Param4 
% External Parameters:
r_tube = Param4(1);
L_tube = Param4(2); %tube length
% T_FG_out6 =Param4 (15);
M_tube = Param4 (5);
n_tube=Param4 (11);
Cp_M = Param4 (16);


% States
Ts_out6= x(1);

% Inputs
RB_p = u(1)/14.7; %bar
steam_flow= u(2); %kg/s from attemporator 2
U = u(3); %220 W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_in6= u(4);
T_FG_out6= u(5);

%properties 
rho_SHRB_in= XSteam ('rho_pT',RB_p,Ts_in6); %density steam (in) from previous compartment
rho_SHRB_out= XSteam ('rho_pT',RB_p, Ts_out6); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',RB_p,Ts_in6) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',RB_p,Ts_out6) ; %Cp steam T_out


%Update Derivative


% Outputs:
Cp_s= (Cp_s_in+Cp_s_out)/2;%Cp steam average
rho_s= (rho_SHRB_in+rho_SHRB_out)/2; %density steam average


% sys(1)is dT_out/dt 

sys(1) =  Ts_out6;

%Add other outputs


% ******************************************************************************************************************************
% Derivatives
% ******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param4)

global  Param4 

% External Parameters:
r_tube = Param4(1);
L_tube = Param4(2); %tube length
% T_FG_out6 =Param4 (14);
M_tube = Param4 (5);
n_tube=Param4 (11);
Cp_M = Param4 (16);


% States
Ts_out6= x(1)

% Inputs
RB_p = u(1)/14.7; %bar
steam_flow= u(2); %kg/s from attemporator 1
U = u(3); %220 W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_in6= u(4);
T_FG_out6= u(5);

%properties 
rho_SHRB_in= XSteam ('rho_pT',RB_p,Ts_in6); %density steam (in) from previous compartment
rho_SHRB_out= XSteam ('rho_pT',RB_p, Ts_out6); %density steam (out) 
Cp_s_in = XSteam ('Cp_pT',RB_p,Ts_in6) ; %Cp steam T_in
Cp_s_out = XSteam ('Cp_pT',RB_p,Ts_out6) ; %Cp steam T_out


% Outputs:
Cp_s= (Cp_s_in+Cp_s_out)/2;%Cp steam average
rho_s= (rho_SHRB_in+rho_SHRB_out)/2; %density steam average
% Derivatives:
% sys(1)is dT_out/dt 

sys(1) =(steam_flow *Cp_s *(Ts_in6-Ts_out6) + U*(2*pi*r_tube*L_tube*n_tube)*(T_FG_out6-Ts_out6))/(rho_s*Cp_s*(pi*r_tube^2*L_tube*n_tube)+M_tube*Cp_M);




