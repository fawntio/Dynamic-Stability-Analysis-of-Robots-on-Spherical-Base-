% clear
% clc

%% ������ �������

drive_params; %�������� ���������� ���������

L = 0.5; %[m]- ����� �������(���������� ����� �������� ����)/����� ��������
g = 10;%[m/c^2] - ��������� ���������� �������
m = 1.5; %[kg]- ����� ����-��������
M = 0.75; %[kg]- ����� ����-���������

R = 0.114;%[m] - ������ ����-��������� (���� 0.12 �)
r = R+0.05;%[m]- ������ ��������

h = abs((L-R)*2);% ������ ��������


%% ��������� �������
h_val = R*2; %[m] - ������ ����-������, ������� ��������, � �������� ������� ���������� �������� R
h_val_korotko = R*2 - R; %[m] - ����� ����������� �����, ��� ��� �� �������� ��� �� ������
r_val = R*0.25; %[m] - ������ ����-������
m_rolik = 1e-5;
%���������� ��� ������������ �������

lenght_rolik = R+r_val; 
static_friction_coeff = 0.05; %����������� ������ �����-������ ���� �� 0.05 �� 0.4
damp_coeff = 1e4;%����������� ��������� ����� ��������
cont_damp = 100;%����������� ������������� �������

% lenght_rolik = sqrt((R + r_val)^2 - R^2);

h = abs((L-R)*2);% ������ ��������

J_rolik = 1/4*m_rolik*r_val^2 + 1/12*m_rolik*h_val^2; %[kg*m^2]- ������ ������� �������� �� ��������� ���� 
J_sphere = 2*M*R^2/5;



%% ���������� ��������� ��������� �������� ������������ � ���������� ��������
syms  x vel tetta tetta_w F Torque xx tetta_ww ang ang_vel;
% syms  L g m M R r J_M J_m;

% x - ��� ���������� ����������� �������
% vel = x'(t) - ��������
% tetta - ��� ���� �������� ��������
% tetta_w = tetta'(t) - ������� �������� ��������


%% �������  �������
J_shar = 2*m*r^2/5; % ������ ������� ����
J_cylinder = 1/4*m*r^2 + 1/12*m*L^2; %[kg*m^2]- ������ ������� �������� �� ��������� ���� 
d = L/2+R;
J_m = J_cylinder + m*d^2;
% ��� d - ���������� ����� ����� ������� (��� ������� ������ �������� �����
%����� ��������)

%������� ������������ � ������������� ������� �������: 
T_base = M / 2 * (vel)^2;% ������������ ������� ���������
T_pend = J_m / 2 * tetta_w^2 + m / 2 * ( (vel)^2 - 2 * L * cos(tetta) * (tetta_w * vel)+...
        L^2*(tetta_w)^2);% ������������ ������� ��������
V_base = 0;% ������������� ������� ���������
V_pend = m*g*L*cos(tetta);% ������������� ������� ��������
%����������
Lag = (T_base+T_pend)-(V_base+V_pend);
% F = m*a*L*cos(tetta+fi); 


%����� ��������� ��������
diff_Lag_vel = M*vel + (m*(2*vel - 2*L*tetta_w*cos(tetta)))/2;
diff_diff_vel_t = M * xx + (m * (2 * xx - 2*L*(tetta_ww*cos(tetta) - tetta_w^2 * sin(tetta)) ) )/2;
diff_Lag_x = 0;

%������ ��������� ��������:
eq_1 = diff_diff_vel_t - diff_Lag_x - Torque/R;


diff_Lag_tetta_w = J_m*tetta_w + (m*(2*tetta_w*L^2 - 2*vel*cos(tetta)*L))/2;
diff_diff_Lag_tetta_w_t = J_m*tetta_ww + (m*(2*tetta_ww*L^2 - 2*L* (xx*cos(tetta) - vel*tetta_w*sin(tetta)) )) / 2;
diff_Lag_tetta = L*g*m*sin(tetta) + L*m*tetta_w*vel*sin(tetta);
%������ ��������� ��������:
eq_2 = diff_diff_Lag_tetta_w_t - diff_Lag_tetta;% ��� Torque - ��� ������, ����������� � ������ ���� ���� (!) �������

%% ��������� ��������� �������
% 
% % �������, ����������� ��� ����������� ��������� ��������� �������
B1 = M + m;
B2 = - m * L * cos(tetta);
B3 = - m * L * cos(tetta);
B4 = J_m + m * L^2;
n1 = m * L * sin(tetta)*tetta_w^2;
n2 = L * tetta_w*vel*sin(tetta) - L * g * m * sin(tetta) - L * m * tetta_w * vel * sin(tetta);

B_q = [B1 B2;B3 B4];
n_q = [n1-Torque/R;n2];

OtVeT = -inv(B_q)*n_q; % ������� ��������� ���������

%%
% ������������

% ��� ������������ ���������� ������� �� ������� ����������� �� ������� ��������� �� ������ �� ����������

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

%������� �� ������� �����������
J = [j1_1 j1_2 j1_3 j1_4;j2_1 j2_2 j2_3 j2_4;j3_1 j3_2 j3_3 j3_4;j4_1 j4_2 j4_3 j4_4];

%%
% ���������� ������� ���������

A = double(subs(J, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));
B = [j1_5; j2_5; j3_5; j4_5];
B = double(subs(B, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));
C = [0 0 1 0];

n = size(A,1);
% �������� �� �������������
r1 = rank(ctrb(A,B));

w = 3; %������ ����������� (�������)
[~,p,~] = besself(r1,w);
% p = [-1; -1.1; -1.2; -1.3]-3.5;


K_reg = place(A,B,p);
Fn = A - B*K_reg;

%% Mass Point

B_q_mass = [M+m -m*L*cos(tetta);-cos(tetta) L];
n_q_mass = [m*L*tetta_w^2*sin(tetta)-Torque/R;-g*sin(tetta)];


OtVeT_mass = -inv(B_q_mass)*n_q_mass; % ������� ��������� ���������

j1_1_mass = diff(vel,x);
j1_2_mass = diff(vel,vel);
j1_3_mass = diff(vel,tetta);
j1_4_mass = diff(vel,tetta_w);
j1_5_mass = diff(vel,Torque);

j2_1_mass = diff((OtVeT_mass(1)),x);
j2_2_mass = diff((OtVeT_mass(1)),vel);
j2_3_mass = diff((OtVeT_mass(1)),tetta);
j2_4_mass = diff((OtVeT_mass(1)),tetta_w);
j2_5_mass = diff((OtVeT_mass(1)),Torque);

j3_1_mass = diff(tetta_w,x);
j3_2_mass = diff(tetta_w,vel);
j3_3_mass = diff(tetta_w,tetta);
j3_4_mass = diff(tetta_w,tetta_w);
j3_5_mass = diff(tetta_w,Torque);

j4_1_mass = diff((OtVeT_mass(2)),x);
j4_2_mass = diff((OtVeT_mass(2)),vel);
j4_3_mass = diff((OtVeT_mass(2)),tetta);
j4_4_mass = diff((OtVeT_mass(2)),tetta_w);
j4_5_mass = diff((OtVeT_mass(2)),Torque);

%������� �� ������� �����������
J_mass = [j1_1_mass j1_2_mass j1_3_mass j1_4_mass;...
    j2_1_mass j2_2_mass j2_3_mass j2_4_mass;...
    j3_1_mass j3_2_mass j3_3_mass j3_4_mass;...
    j4_1_mass j4_2_mass j4_3_mass j4_4_mass];

A_mass = double(subs(J_mass, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));

B_mass = [j1_5_mass; j2_5_mass; j3_5_mass; j4_5_mass];
B_mass = double(subs(B_mass, {x, vel, tetta, tetta_w, Torque}, {0,0,0,0,0}));
C_mass = [1 0 0 0];

r2 = rank(ctrb(A_mass,B_mass));

n_mass = size(A_mass,1);
[~,p_mass,~] = besself(r2,w);
K_reg_mass = place(A_mass,B_mass,p_mass);
Fn_mass = A_mass - B_mass*K_reg_mass;




