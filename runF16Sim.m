%================================================
%     Matlab Script File used to run the
%     non-linear F-16 Simulation.  The results
%     will also be saved to a file and plotted.
%
% Author: Richard S. Russell
% Modified by Wenchao Lei at 20180714
%
%================================================

%% Start Timer
%%
tic;

%% Prepare
%%
clc; clear; close all;

global fi_type fi_flag_Simulink;
global altitude velocity;
global surface1 surface2 surface3;
global ElevatorDis AileronDis RudderDis;

surface1 = 'ele_';
surface2 = 'ail_';
surface3 = 'rud_';

disp('This is an F-16 Simulation.');
disp('The simulation will begin by asking you for the flight ');
disp('conditions for which the simulation will be performed.');
disp(newline);
disp('Accpetable values for flight condition parameters are:');
disp(newline);
disp('                                  Model');
disp('  Variable                LOFI            HIFI');
disp('              Units   Min     Max     Min     Max');
disp('  Altitude:   ft      5000    40000   5000    40000');
disp('  AOA         deg    -10      45     -10      90');
disp('  Thrust      lbs     1000    19000   1000    19000');
disp('  Elevator    deg    -25.0    25.0   -25.0    25.0');
disp('  Aileron     deg    -21.5    21.5   -21.5    21.5');
disp('  Rudder      deg    -30      30     -30      30');
disp('  Velocity    ft/s    300     900     300     900');
disp(newline);
disp('The flight condition you choose will be used to trim the F16.');
disp('Note:  The trim routine will trim to the desired');
disp('altitude and velocity.  All other parameters');
disp('will be varied until level flight is achieved.  ');
disp('You may need to view the results of the simulation');
disp(' and retrim accordingly.');
disp(newline);

%% Ask user which simulation to run.
%%
disp('Which model would you like to use to trim the aircraft:')
disp('  1. Low Fidelity F-16 Trim')
disp('  2. High Fidelity F-16 Trim')
fi_flag = input('Your Selection:  ');
disp(newline);

%% Determine from flag the correct simulation.
%%
if fi_flag == 1
    fi_type = 'lofi';
    fi_flag_Simulink = 0;
elseif fi_flag == 2
    fi_type = 'hifi';
    fi_flag_Simulink = 1;
else
    disp('Invalid selection and we will choose the default mode: HIFI Model');
    disp(newline);
    fi_type = 'hifi';
    fi_flag_Simulink = 1;
end

%% Trim aircraft to desired altitude and velocity
%%
altitude = input('Enter the altitude for the simulation (ft)  :  ');
velocity = input('Enter the velocity for the simulation (ft/s):  ');

%% Initialize some varibles used to create disturbances.
%%
DisEle_1 = 0;    DisEle_2 = 0;    DisEle_3 = 0;
DisAil_1 = 0;    DisAil_2 = 0;    DisAil_3 = 0;
DisRud_1 = 0;    DisRud_2 = 0;    DisRud_3 = 0;
ElevatorDis = 0; AileronDis = 0;  RudderDis = 0;

%% Find out which surface to create a disturbance on.
%%
dis_flag = input('Would you like to create a disturbance on a surface (y/n):  ', 's');
disp(newline);

if dis_flag == 'y'
    ElevatorDis = input('Enter the elevator distrubance deflection    (deg) :  ');
    DisEle_1 = ElevatorDis;
    DisEle_2 = -2*ElevatorDis;
    DisEle_3 = ElevatorDis;
    
    AileronDis = input('Enter the aileron distrubance deflection      (deg) :  ');
    DisAil_1 = AileronDis;
    DisAil_2 = -2*AileronDis;
    DisAil_3 = AileronDis;
    
    RudderDis = input('Enter the rudder distrubance deflection        (deg) :  ');
    DisRud_1 = RudderDis;
    DisRud_2 = -2*RudderDis;
    DisRud_3 = RudderDis;
    
    surfacedis = 'ele_ail_rud';
    
elseif dis_flag == 'n'
    surfacedis = 'none';
    
else
    disp('Invalid Selection and we will choose the default value: n!');
    disp(newline);
    surfacedis = 'none';
end

%% Time Setting
%%
delta_T = 0.001;
TStart = 0; TFinal = 60;

%% Initial Conditions for trim routine.
%================================================
% The following values seem to trim to most flight condition.
% If the F16 does not trim, change these values.
%================================================
%%
alpha = 8.49;           % AOA, degrees
thrust = 5000;          % thrust, lbs
elevator = -0.09;       % elevator, degrees
rudder = -0.01;         % rudder angle, degrees
aileron = 0.01;         % aileron, degrees
%% Trim the Initial Conditions
%%
[trim_state, trim_thrust, trim_control, dLEF, UX] = trim_F16(altitude, velocity, alpha, thrust, elevator, aileron, rudder);

%% Simulate the F-16 Dynamics
%%
%trim_control = [0;0;0];
sim( 'F16Block' ,[TStart TFinal]);

%% Prepare the file to save the data
%%
trim_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, fi_type, altitude, velocity);
fid_trim = fopen(trim_file, 'w');
%heading1 = sprintf('\n\t\t  %s DATA Trim-Doublet on %s: Alt %.0f, Alpha %.0f\n\n', fi_type, surfacedef, altitude, alpha);
heading2 = sprintf('\ntime,npos,epos,alt,phi,theta,psi,vel,alpha,beta,p,q,r,nx,ny,nz,mach,qbar,ps,\n\n');
%fprintf(fid_trim,heading1);
fprintf(fid_trim,heading2);

%% Save the simulation data
%%
fid_trim = fopen(trim_file, 'a');
for row = 1 : 1 : length(y_sim(:,1))
    fprintf(fid_trim,'%8.5f,',T(row,:));
    for column = 1 : 1 : length(y_sim(1,:))
        fprintf(fid_trim,'%8.5f,',y_sim(row,column));
    end
    for column = 1:1:length(surfaces(1,:))
        fprintf(fid_trim,'%8.5f,',surfaces(row,column));
    end
    fprintf(fid_trim,'\n');
end
fclose(fid_trim);

%% Plot
%%
plot_flag = input('Plot results (y/n):  ', 's');
if plot_flag == 'y'
    graphF16;
end

%% End Timer
%%
toc
