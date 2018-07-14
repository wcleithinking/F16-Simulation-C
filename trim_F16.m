function [trim_state, trim_thrust, trim_control, dLEF, xu] = trim_F16(alt, vel, alpha, thrust, elevator, ail, rud)
%================================================
%     F16 nonlinear model trimming routine
%  for longitudinal motion, steady level flight
%
% Author: T. Keviczky
% Date:   April 29, 2002
%
%
%      Added addtional functionality.
%      This trim function can now trim at three
%      additional flight conditions
%         -  Steady Turning Flight given turn rate
%         -  Steady Pull-up flight - given pull-up rate
%         -  Steady Roll - given roll rate
%
% Coauthor: Richard S. Russell
% Date:     November 7th, 2002
%
% Modified by Wenchao Lei at 20180714
%
%================================================
% OUTPUTS: trimmed values for states and controls
% INPUTS:  guess values for thrust, elevator, alpha  (assuming steady level flight)
%================================================

% global fi_flag_Simulink;
global altitude velocity;
global phi psi;
global p q r;
global phi_weight theta_weight psi_weight

altitude = alt;
velocity = vel;
alpha = alpha*pi/180;  % convert to radians

% Initial Guess for free parameters
UX0 = [alpha; thrust; elevator; ail; rud];  % free parameters: two control values & angle of attack

% Initialize some varibles
phi = 0; psi = 0;
p = 0; q = 0; r = 0;
phi_weight = 10; theta_weight = 10; psi_weight = 10;

% modify the variables according to the given flight condition
disp('At what flight condition would you like to trim the F-16?');
disp('1.  Steady Wings-Level Flight.');
disp('2.  Steady Turning Flight.');
disp('3.  Steady Pull-Up Flight.');
disp('4.  Steady Roll Flight.');
FC_flag = input('Your Selection:  ');

switch FC_flag
    case 1
        % do nothing
    case 2
        r = input('Enter the turning rate (deg/s):  ');
        psi_weight = 0;
    case 3
        q = input('Enter the pull-up rate (deg/s):  ');
        theta_weight = 0;
    case 4
        p = input('Enter the Roll rate    (deg/s):  ');
        phi_weight = 0;
    otherwise
        disp('Invalid Selection and we will choose the default value: <1. Steady Wings-Level Flight.>!!!')
end

% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+04,'MaxIter',1e+04);

iter = 1;
while iter == 1
    
    %[UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('trimfun',UX0,OPTIONS);
    UX = fminsearch('trimfun',UX0,OPTIONS);
    
    [cost, ~, xu] = trimfun(UX);
    
    disp('Trim Values and Cost:');
    disp(['cost   = ' num2str(cost)])
    disp(['thrust = ' num2str(xu(13)) ' lb'])
    disp(['elev   = ' num2str(xu(14)) ' deg'])
    disp(['ail    = ' num2str(xu(15)) ' deg'])
    disp(['rud    = ' num2str(xu(16)) ' deg'])
    disp(['alpha  = ' num2str(xu(8)*180/pi) ' deg'])
    disp(['dLEF   = ' num2str(xu(17)) ' deg'])
    disp(['Vel.   = ' num2str(velocity) 'ft/s'])
    flag = input('Continue trim rountine iterations? (y/n):  ','s');
    if flag == 'n'
        iter = 0;
    end
    UX0 = UX;
end

% For simulink:
trim_state=xu(1:12);
trim_thrust=UX(2);
trim_control=[UX(3);UX(4);UX(5)];
dLEF = xu(17);