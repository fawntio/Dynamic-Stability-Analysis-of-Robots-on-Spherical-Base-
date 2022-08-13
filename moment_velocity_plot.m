%% Plotting
%Settings
set(0,'defaultAxesFontSize',13);
set(0,'defaultTextFontSize',13);
set(groot, 'defaultFigureUnits', 'normalized');
set(groot, 'defaultFigurePosition', [0.2 0.2 0.5 0.6]);
set(groot, 'defaultLegendLocation', 'best');
set(groot,'defaultAxesNextPlot','add');
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 

%% Running the Model
Iitial_code; %running the code with the calculation of the coefficients of the regulator 
Torque = sim('roliki_model'); % running the model 

%% Plotting
% Moment/Torque
fig = figure('Name','Graph of moments on the roller and on the ball');
plot(torque_x.time, torque_x.data,'r-','DisplayName', 'f = M_{sphere}')
hold on
plot(torque_x.time,torque_x.data*(-r_val/R) ,'c--','DisplayName', 'f = M_{rolik}')
grid on
legend
yline([0.03])
yline([-0.03])

xlabel('t, s');
ylabel('Torque_x, N*m');
saveas(fig, strcat(get(fig, 'Name'), '.png'));
saveas(fig, strcat(get(fig, 'Name'), '.fig'));
close(fig);
%%
% Velocity
fig = figure('Name','Graph of the speed on the roller and on the ball');
plot(sp_rad_vel_x.time, sp_rad_vel_x.data,'r-','DisplayName', 'f = v_{sphere}')
hold on
plot(sp_rad_vel_x.time,sp_rad_vel_x.data*(-R/r_val) ,'c--','DisplayName', 'f = v_{rolik}')
grid on
legend
xlabel('t, s');
ylabel('Vel_x, rad/s');
saveas(fig, strcat(get(fig, 'Name'), '.png'));
saveas(fig, strcat(get(fig, 'Name'), '.fig'));
close(fig);

fig = figure('Name','Graph of the speed on the roller and on the ball in rpm');
plot(sp_rad_vel_x.time, sp_rad_vel_x.data*2*pi,'r-','DisplayName', 'f = v_{sphere}')
hold on
plot(sp_rad_vel_x.time,sp_rad_vel_x.data*(-R/r_val)*2*pi,'c--','DisplayName', 'f = v_{rolik}')
grid on
legend
xlabel('t, s');
ylabel('Vel_x, rev per minute');
saveas(fig, strcat(get(fig, 'Name'), '.png'));
saveas(fig, strcat(get(fig, 'Name'), '.fig'));
close(fig);

