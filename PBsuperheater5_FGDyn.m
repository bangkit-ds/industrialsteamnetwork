function [sys,x0,str,tss]=PBsuperheater5_FGDyn(t,x,u,flag,Param5,X_ss)

global  Param5
switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
        % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param5);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param5);

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

global  Param5 

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


function sys = mdlOutputs(t,x,u,Param5);

global  Param5
%External Parameters:
r_tube = Param5(1);
L_tube = Param5(2); %tube length
M_tube = Param5 (3);
n_tube=Param5 (4);
Cp_FG = Param5 (5);
Cp_M = Param5 (6);
A_FG = Param5 (7);
rho_FG= Param5 (8);


% States
TFG_out5= x(1);

% Inputs
FG_flow= u(1)/3600; %kg/s
T_FGin5= u(2); %oC
U = u(3); % W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_out5= u(4); %oC





% Outputs Auxilary:
Q_5=FG_flow*Cp_FG *(T_FGin5-TFG_out5);

%Output
% states

sys(1) =  TFG_out5;

%Auxiliary variables
sys(2)= Q_5;


% ******************************************************************************************************************************
% Derivatives
%%******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param5)

global  Param5 

%External Parameters:
r_tube = Param5(1);
L_tube = Param5(2); %tube length
M_tube = Param5 (3);
n_tube=Param5 (4);
Cp_FG = Param5 (5);
Cp_M = Param5 (6);
A_FG = Param5 (7);
rho_FG= Param5 (8);



%States
TFG_out5= x(1);

%Inputs
FG_flow= u(1)/3600; %kg/s
T_FGin5= u(2); %oC
U = u(3); % W/m2K -> make it input and simulate for fouling make it ramping and became prediction for soot blowing activities
Ts_out5= u(4); %oC



% Outputs Auxilary:
Q_5=FG_flow*Cp_FG *(T_FGin5-TFG_out5);
% 
% 
% Derivatives: TFG_out6

sys(1) =(FG_flow*(Cp_FG* T_FGin5 - Cp_FG*TFG_out5) - U*(2*pi*r_tube*L_tube*n_tube)*(TFG_out5-Ts_out5))/(rho_FG*Cp_FG*(A_FG*L_tube*n_tube)+M_tube *Cp_M);



