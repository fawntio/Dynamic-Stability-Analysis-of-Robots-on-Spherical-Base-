
%% Estimated real parameters
L = 0.3; %[m]- suspension length (distance between centers of mass)
g = 10;%[m/c^2] - acceleration of gravity
m = 1.5; %[kg]- pendulum mass
M = 0.75; %[kg]- spherical base mass

R = 0.114;%[m] - spherical base radius
r = R+0.05;%[m]- pendulum radius

%% Параметры роликов
h_val_korotko = R*2 - R; %[m] - roller length
r_val = R*0.25; %[m] - roller radius
m_rolik = 1e-5;%[kg]- roller mass

lenght_rolik = R+r_val; %[m] - roller spacing
static_friction_coeff = 0.05; %friction coefficient steel-rubber goes from 0.05 to 0.4
damp_coeff = 1e4;%coefficient of elasticity between rollers
cont_damp = 100;%roller damping factor

h = abs((L-R)*2);% pendulum length

J_rolik = 1/4*m_rolik*r_val^2 + 1/12*m_rolik*h_val^2; %[kg*m^2]- moment of inertia of an offset cylinder
J_sphere = 2*M*R^2/5;


%% We compose the equations of system behavior in accordance with the Lagrange equation
syms  x vel tetta tetta_w F Torque xx tetta_ww ang ang_vel al Om V;
% syms  L g m M R r J_M J_m;

% x - coordinate of the carriage movement
% vel = x'(t) - speed
% tetta - pendulum rotation angle
% tetta_w = tetta'(t) - pendulum angular velocity
%% Moments of inertia
J_shar = 2*m*r^2/5; % moment of inertia of the ball
J_cylinder = 1/4*m*r^2 + 1/12*m*h^2; %[kg*m^2]- moment of inertia of an offset cylinder
d = L/2-R;
% d = L;
J_m = J_cylinder + m*d^2;
% where d is the distance between the axes of inertia (the second axis of inertia passes through
% end of cylinder)

%% Equations

%We express the kinetic and potential energies of the system:
T_base = M / 2 * (vel)^2;% base kinetic energy
T_pend = J_m / 2 * tetta_w^2 + m / 2 * ( (vel)^2 - 2 * L * cos(tetta) * (tetta_w * vel)+...
        L^2*(tetta_w)^2);% pendulum kinetic energy
V_base = 0;% base potential energy
V_pend = m*g*L*cos(tetta);% potential energy of the pendulum
%Lagrangian
Lag = (T_base+T_pend)-(V_base+V_pend);

F = Torque/R;
%Derivation of the equations of motion
diff_Lag_vel = M*vel + (m*(2*vel - 2*L*tetta_w*cos(tetta)))/2;
diff_diff_vel_t = M * xx + (m * (2 * xx - 2*L*(tetta_ww*cos(tetta) - tetta_w^2 * sin(tetta)) ) )/2;
diff_Lag_x = 0;

%First equation of motion:
eq_1 = diff_diff_vel_t - diff_Lag_x - Torque/R;


diff_Lag_tetta_w = J_m*tetta_w + (m*(2*tetta_w*L^2 - 2*vel*cos(tetta)*L))/2;
diff_diff_Lag_tetta_w_t = J_m*tetta_ww + (m*(2*tetta_ww*L^2 - 2*L* (xx*cos(tetta) - vel*tetta_w*sin(tetta)) )) / 2;
diff_Lag_tetta = L*g*m*sin(tetta) + L*m*tetta_w*vel*sin(tetta);
%Second equation of motion:
eq_2 = diff_diff_Lag_tetta_w_t - diff_Lag_tetta;% где Torque - это момент, приложенный к центру масс всей (!) системы

%% System State Equations
% Matrices necessary for compiling the equations of state of the system
B1 = M + m;
B2 = - m * L * cos(tetta);
B3 = - m * L * cos(tetta);
B4 = J_m + m * L^2;
n1 = m * L * sin(tetta)*tetta_w^2;
n2 = L * tetta_w*vel*sin(tetta) - L * g * m * sin(tetta) - L * m * tetta_w * vel * sin(tetta);

B_q = [B1 B2;B3 B4];
n_q = [n1-Torque/R;n2];

OtVeT = -inv(B_q)*n_q; % state equations matrix

%%
% Linearization

% For linearization, we compose a matrix of partial derivatives for each equation for each of the variables
j1_1 = diff(vel,x);
j1_2 = diff(vel,vel);
j1_3 = diff(vel,tetta);
j1_4 = diff(vel,tetta_w);
j1_5 = diff(vel,Torque);

j2_1 = diff((OtVeT(1)),x);
j2_2 = diff((OtVeT(1)),vel);
j2_3 = diff((OtVeT(1)),tetta);
j2_4 = diff((OtVeT(1)),tetta_w);
j2_5 = diff((OtVeT(1)),Torque);

j3_1 = diff(tetta_w,x);
j3_2 = diff(tetta_w,vel);
j3_3 = diff(tetta_w,tetta);
j3_4 = diff(tetta_w,tetta_w);
j3_5 = diff(tetta_w,Torque);

j4_1 = diff((OtVeT(2)),x);
j4_2 = diff((OtVeT(2)),vel);
j4_3 = diff((OtVeT(2)),tetta);
j4_4 = diff((OtVeT(2)),tetta_w);
j4_5 = diff((OtVeT(2)),Torque);

%Matrix of partial derivatives
J = [j1_1 j1_2 j1_3 j1_4;j2_1 j2_2 j2_3 j2_4;j3_1 j3_2 j3_3 j3_4;j4_1 j4_2 j4_3 j4_4];

%%
% Compiling state matrices

A = double(subs(J, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));
B = [j1_5; j2_5; j3_5; j4_5];
B = double(subs(B, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));
C = [1 0 0 0];

n = size(A,1);
% Controllability test
r1 = rank(ctrb(A,B));

w = 3; %bandwidth (frequency)
[~,p,~] = besself(r1,w); %bessel filter
p = [-1; -1.1; -1.2; -1.3]-3.5;


K_reg = place(A,B,p);
Fn = A - B*K_reg;

%% 
% Creating a regulator with an internal system

Gr = [0 1;0 0]; 
Hr = [1 0]; 
Lr = [0; 1];
nr = size(Gr,1); 

% State matrices for a controller with an internal system
An = [A zeros(n, nr);-Lr*C Gr];
Bn = [B; zeros(nr,1)];
Cn = [C zeros(1,nr)];

[~,p,~] = besself(n+nr, w);
% p = [-1; -1.1; -1.2; -1.3; -1.4; -1.5]-3.5;

K = place(An,Bn,p);
Fn = An - Bn*K; % closed system matrix

% Controllability test
r2 = rank(ctrb(An,Bn));

eig(Fn);



