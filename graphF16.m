function graphF16()

%% Function to Graph the results of runsim
%%

global fi_type;
global altitude velocity;
global surface1 surface2 surface3;
global ElevatorDis AileronDis RudderDis;

lofi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'lofi', altitude, velocity);
hifi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'hifi', altitude, velocity);

lofiID = fopen(lofi_data_file,'r');
hifiID = fopen(hifi_data_file,'r');

title_string = sprintf('Trimed at Velocity = %.1f \n Alt. = %.1f', velocity, altitude);

if (lofiID > 0 && hifiID > 0)
    
    C_lofi = textscan(lofiID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    C_hifi = textscan(hifiID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    fclose(lofiID);
    fclose(hifiID);
    time_lo  = C_lofi{1};
    npos_lo  = C_lofi{2};
    epos_lo  = C_lofi{3};
    alt_lo   = C_lofi{4};
    phi_lo   = C_lofi{5};
    theta_lo = C_lofi{6};
    psi_lo   = C_lofi{7};
    vel_lo   = C_lofi{8};
    alpha_lo = C_lofi{9};
    sideslip_lo = C_lofi{10};
    roll_lo  = C_lofi{11};
    pitch_lo = C_lofi{12};
    yaw_lo   = C_lofi{13};
    nx_lo    = C_lofi{14};
    ny_lo    = C_lofi{15};
    nz_lo    = C_lofi{16};
    mach_lo  = C_lofi{17};
    qbar_lo  = C_lofi{18};
    ps_lo    = C_lofi{19};
    thrust_lo = C_lofi{20};
    ele_lo   = C_lofi{21};
    ail_lo   = C_lofi{22};
    rud_lo   = C_lofi{23};
    
    time_hi  = C_hifi{1};
    npos_hi  = C_hifi{2};
    epos_hi  = C_hifi{3};
    alt_hi   = C_hifi{4};
    phi_hi   = C_hifi{5};
    theta_hi = C_hifi{6};
    psi_hi   = C_hifi{7};
    vel_hi   = C_hifi{8};
    alpha_hi = C_hifi{9};
    sideslip_hi = C_hifi{10};
    roll_hi  = C_hifi{11};
    pitch_hi = C_hifi{12};
    yaw_hi   = C_hifi{13};
    nx_hi    = C_hifi{14};
    ny_hi    = C_hifi{15};
    nz_hi    = C_hifi{16};
    mach_hi  = C_hifi{17};
    qbar_hi  = C_hifi{18};
    ps_hi    = C_hifi{19};
    thrust_hi = C_hifi{20};
    ele_hi   = C_hifi{21};
    ail_hi   = C_hifi{22};
    rud_hi   = C_hifi{23};
    
    %[time_lo, npos_lo, epos_lo, alt_lo, phi_lo, theta_lo, psi_lo, vel_lo, alpha_lo, sideslip_lo, roll_lo, pitch_lo, yaw_lo, nx_lo, ny_lo, nz_lo, mach_lo, qbar_lo, ps_lo, thrust_lo, ele_lo, ail_lo, rud_lo] = textscan(lofi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    %[time_hi, npos_hi, epos_hi, alt_hi, phi_hi, theta_hi, psi_hi, vel_hi, alpha_hi, sideslip_hi, roll_hi, pitch_hi, yaw_hi, nx_hi, ny_hi, nz_hi, mach_hi, qbar_hi, ps_hi, thrust_hi, ele_hi, ail_hi, rud_hi] = textscan(hifi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    
    %% Figure 1
    %%
    figure(1);
    subplot(231)
    plot(time_lo,npos_lo, time_hi, npos_hi, '--');
    ylabel('North Pos.');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(232)
    plot(time_lo, epos_lo, time_hi, epos_hi, '--');
    ylabel('East Pos.');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(233)
    plot(time_lo, alt_lo,time_hi , alt_hi, '--') ;
    ylabel('Altitude ft');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(234)
    plot(time_lo, phi_lo, time_hi, phi_hi, '--');
    ylabel('PHI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(235)
    plot(time_lo, theta_lo, time_hi, theta_hi, '--');
    ylabel('THETA (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(236)
    plot(time_lo ,psi_lo, time_hi, psi_hi, '--');
    ylabel('PSI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 2
    %%
    
    figure(2)
    subplot(231)
    plot(time_lo,vel_lo, time_hi, vel_hi, '--');
    ylabel('Velocity (ft/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(232)
    plot(time_lo, alpha_lo, time_hi, alpha_hi, '--');
    ylabel('Angle of Attack (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(233)
    plot(time_lo , sideslip_lo, time_hi, sideslip_hi, '--');
    ylabel('Side Slip (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(234)
    plot(time_lo,roll_lo, time_hi, roll_hi, '--');
    ylabel('Roll Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(235)
    plot(time_lo, pitch_lo , time_hi, pitch_hi, '--');
    ylabel('Pitch Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(236)
    plot(time_lo, yaw_lo, time_hi, yaw_hi, '--');
    ylabel('Yaw Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 3
    %%
    figure(3);
    subplot(231)
    plot(time_lo, nx_lo, time_hi, nx_hi, '--');
    ylabel('acc x');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(232)
    plot(time_lo, ny_lo, time_hi, ny_hi, '--');
    ylabel('acc y');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(233)
    plot(time_lo , nz_lo, time_hi, nz_hi, '--');
    ylabel('acc z');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(234)
    plot(time_lo, mach_lo, time_hi, mach_hi, '--');
    ylabel('Mach');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(235)
    plot(time_lo, qbar_lo, time_hi, qbar_hi, '--');
    ylabel('q bar)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(236)
    plot(time_lo, ps_lo, time_hi, ps_hi, '--');
    ylabel('ps');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 4
    %%
    figure(4)
    subplot(221)
    plot(time_lo, thrust_lo, time_hi, thrust_hi, '--');
    ylabel('del Thrust');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(222)
    plot(time_lo, ele_lo, time_hi, ele_hi, '--');
    ylabel('del Elevator');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(223)
    plot(time_lo , ail_lo, time_hi, ail_hi, '--');
    ylabel('del Aileorn');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(224)
    plot(time_lo, rud_lo, time_hi, rud_hi,'--');
    ylabel('del Rudder');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 5
    %%
    alt_min_hi = min(alt_hi) - 2000;
    alt_min(1:length(alt_lo),1) = alt_min_hi;
    figure(5)
    plot3(epos_lo,npos_lo,alt_lo, 'b-' , epos_hi,npos_hi,alt_hi,'g-.', epos_lo,npos_lo,alt_min, 'b:',  epos_hi,npos_hi,alt_min,'g:','linewidth',2)
    grid;
    legend('LOFI', 'HIFI','LOFI_{xy}','HIFI_{xy}');
    xlabel('East Position')
    ylabel('North Position')
    zlabel('Altitude')
    
    %% Figure 6
    %%
    figure(6)
    for j = 1:10:length(epos_lo)
        plot3(epos_lo,npos_lo,alt_lo, 'b:',epos_hi,npos_hi,alt_hi, 'g:', epos_lo(j),npos_lo(j),alt_lo(j), 'ro', epos_hi(j),npos_hi(j),alt_hi(j), 'ro','linewidth',2)
        grid;
        xlabel('East Position')
        ylabel('North Position')
        zlabel('Altitude')
        legend('LOFI', 'HIFI')
        F(j) = getframe;
    end
    
else
    
    new_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, fi_type, altitude, velocity);
    newID = fopen(new_data_file,'r');
    C_new = textscan(newID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    fclose(newID);
    time_new  = C_new{1};
    npos_new  = C_new{2};
    epos_new  = C_new{3};
    alt_new   = C_new{4};
    phi_new   = C_new{5};
    theta_new = C_new{6};
    psi_new   = C_new{7};
    vel_new   = C_new{8};
    alpha_new = C_new{9};
    sideslip_new = C_new{10};
    roll_new  = C_new{11};
    pitch_new = C_new{12};
    yaw_new   = C_new{13};
    nx_new    = C_new{14};
    ny_new    = C_new{15};
    nz_new    = C_new{16};
    mach_new  = C_new{17};
    qbar_new  = C_new{18};
    ps_new    = C_new{19};
    thrust_new = C_new{20};
    ele_new   = C_new{21};
    ail_new   = C_new{22};
    rud_new   = C_new{23};
    %[time_new, npos_new, epos_new, alt_new, phi_new, theta_new, psi_new, vel_new, alpha_new, sideslip_new, roll_new, pitch_new, yaw_new, nx_new, ny_new, nz_new, mach_new, qbar_new, ps_new, thrust_new, ele_new, ail_new, rud_new] = textread(new_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    
    %% Figure 1
    %%
    figure(1);
    subplot(231)
    plot(time_new,npos_new);
    ylabel('North Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, epos_new);
    ylabel('East Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , alt_new);
    ylabel('Altitude ft');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new, phi_new);
    ylabel('Roll: \phi (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new, theta_new);
    ylabel('Pitch: \theta (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new ,psi_new);
    ylabel('Yaw: \psi (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 2
    %%
    figure(2);
    subplot(231)
    plot(time_new,vel_new);
    ylabel('Velocity (ft/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, alpha_new);
    ylabel('Angle of Attack: \alpha (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , sideslip_new);
    ylabel('Side Slip: \beta (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new,roll_new);
    ylabel('Roll Rate: p (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new,pitch_new);
    ylabel('Pitch Rate: q (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new, yaw_new);
    ylabel('Yaw Rate: r (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 3
    %%
    figure(3);
    subplot(231)
    plot(time_new, nx_new);
    ylabel('acc x');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, ny_new);
    ylabel('acc y');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , nz_new);
    ylabel('acc z');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new, mach_new);
    ylabel('Mach');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new, qbar_new);
    ylabel('q bar)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new, ps_new);
    ylabel('ps');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 4
    %%
    figure(4);
    subplot(221)
    plot(time_new, thrust_new);
    ylabel('Thrust');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(222)
    plot(time_new, ele_new);
    ylabel('del Elevator');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(223)
    plot(time_new , ail_new);
    ylabel('del Aileron');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(224)
    plot(time_new, rud_new);
    ylabel('del Rudder');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 5
    %%
    alt_min_new = min(alt_new) - 2000;
    alt_min(1:length(alt_new),1) = alt_min_new;
    figure(5)
    plot3(epos_new, npos_new, alt_new, 'b',epos_new,npos_new,alt_min, 'b:','linewidth',2)
    grid;
    if strcmpi(fi_type,'lofi')
        legend('LOFI', 'LOFI_{xy}');
    elseif strcmpi(fi_type,'hifi')
        legend('HIFI', 'HIFI_{xy}');
    end
    xlabel('East Position')
    ylabel('North Position')
    zlabel('Altitude')
    
    %% Figure 6
    %%
    figure(6)
    for j = 1:10:length(epos_new)
        plot3(epos_new, npos_new, alt_new, 'b:', epos_new(j),npos_new(j),alt_new(j), 'ro','linewidth',2)
        grid;
        xlabel('East Position')
        ylabel('North Position')
        zlabel('Altitude')
        if strcmpi(fi_type,'lofi')
            legend('LOFI');
        elseif strcmpi(fi_type,'hifi')
            legend('HIFI');
        end
        F(j) = getframe;
    end
end  % end if

movie(F)

end