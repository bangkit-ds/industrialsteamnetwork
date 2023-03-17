function [sys,x0,str,tss]=cstr1(t,x,u,flag,Param,X_ss)

global  ph



switch flag,

case 0,	% Initialize the states and sizes
   [sys,x0,str,tss] = mdlInitialSizes(t,x,u,X_ss);
   
        % ****************
  	%  Outputs
  	% ****************   
     
case 3,   % Calculate the outputs
   
   sys = mdlOutputs(t,x,u,Param);
   
   % ****************
   % Update
   % ****************
case 1,	% Obtain derivatives of states
   
   sys = mdlDerivatives(t,x,u,Param);

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

global  ph

% This handles initialization of the function.
% Call simsize of a sizes structure.
sizes = simsizes;
sizes.NumContStates  = 1;     % continuous states
sizes.NumDiscStates  = 0;     % discrete states
sizes.NumOutputs     = 1;    % outputs
sizes.NumInputs      = 4;     % inputs 
sizes.DirFeedthrough = 0;     % System is non-causal
sizes.NumSampleTimes = 1;     %
sys = simsizes(sizes);        %
x0  = X_ss;                   % Initialize the continuous states.
str = [];                     % set str to an empty matrix.
tss = [0,0];	              % sample time: [period, offset].















% ******************************************************************************************************************************
%  Outputs
% ******************************************************************************************************************************


function sys = mdlOutputs(t,x,u,Param);

global  ph

% External Parameters:

%Thdr    = Param(1);
%Vhdr    = Param(2);


% States

ph = x(1);


% Inputs

dqRB = u(1);
dqPB = u(2);
dqSU = u(3);
dqTU = u(4);



% Outputs:

% States

sys(1) = ph;



% ******************************************************************************************************************************
% Derivatives
% ******************************************************************************************************************************

function sys = mdlDerivatives(t,x,u,Param)

global ph

% External Parameters:

%Thdr    = Param(1);
%Vhdr    = Param(2);


% States

ph  = x(1);

% Inputs

dqRB = u(1);
dqPB = u(2);
dqSU = u(3);
dqTU = u(4);




% Derivatives:


% sys(1)is dph/dt 

sys(1) = 5.581157407407407e+004/3*1e-6*(dqRB+dqPB-dqSU-dqTU);
