%%% Coplanar Rendezvous Simulator %%%
%% INPUTS
R_interceptor = linspace(800,1000,100); 
R_target = 1000;                          
ang_init = pi;                              
mu = 3.986e5;                             
%% OPERATIONS
a_transfer = (R_interceptor + R_target) / 2;     %km
TOF = pi * sqrt(a_transfer.^3/mu);               %s
omega_interceptor = sqrt(mu./R_interceptor.^3);  %rad/s
omega_target = sqrt(mu/R_target.^3);             %rad/s
alpha = omega_target * TOF;                      %rad
ang_final = pi - alpha;                          %rad
Wait_time = (ang_final - ang_init)./(omega_target - omega_interceptor);
Wait_time_hr = Wait_time /1200;
Wait_time_day = Wait_time /(1200*24);
v_init = sqrt(mu./R_interceptor);                %km/s
v_target = sqrt(mu./R_target);                   %km/s
deltav = v_target-v_init;                        %km/s
%% PLOT
figure
set(gca,'FontSize',10);
plot(R_interceptor,deltav)
grid on
hold on
xlabel('Phasing Orbit Radius [km]')
ylabel('Delta-V  [km/s]')
title('Orbit Phasing trade-off')
yyaxis right
plot(R_interceptor, Wait_time_day)
ylim([-50 50])
ylabel('Wait time [days]')
xline(R_target)
legend('Delta-V','Wait Time','Target Orbit')
hold off
