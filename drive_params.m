%% Параметры двигателя

L_dr = 0.5e-3; 
R_dr = 5.5; 
U_nom = 12; %[В] - номинальное напряжение
Udc = U_nom*sqrt(2);
Js = 1.0519*10^(-7);
n=2375;
Ce=0.0527;
% Ce=0.3098;
Cm = Ce;
% Cm=0.0028;

% 
% Te = L_dr/R_dr;
% Tt = 0.001;
% Ki = R_dr/Tt;
% Kp = Te*Ki;

wn=2*pi*n/60;
Cea=Ce/sqrt(2);  %В*с/рад, постоянная противоЭДС двигателя А
Cma=Cea;         %Нм/А, коэффициент момента двигателя А

Ra=R_dr/2;         %Активное сопротивление фазы двигателя А
La=L_dr/2;         %Индуктивность фазы двигателя А

Te=1.5*La/Ra;
Ld=La;
Lq=Ld;

Isc=Udc/(Ra*3/2); % Ток КЗ
Msc=Cma*Isc;     % Максимальный момент


%% Параметры преобразователя

%Частота коммутации fk
Ts=5e-6;
fk=100000 ;                       %   -   Частота коммутации [Гц]
Tk=1/fk ;                         %   -   Период коммутации [с]
T=Tk;

% Tconv=T;                % двухсторонняя пила 0.5 % - Converter time constant, s;
Tconv=0.5*T;
% Kconv=1/wn;
Kconv=Udc;                 % - Converter koefficient ;
b=Cma*3/2*Cea/Ra;        % - Speed-Torque characteristic stiffness;
% b=Cm*I_Md/dw
% НАСТРОЙКА КОНТУРА ТОКА (технический оптимум)
In=0.17*sqrt(2);
Ti=Te;                                   % - Current regulator time constant
Kfb_i=1/In;                              % - Current feedback coefficient
Ki=Te/(2*Kfb_i*(1/Cma)*Tconv*Kconv/Cea*b);   % - current regulator gain
% Ki=Te/(2*T*Kconv*b)                   % - current regulator gain in pu
Rshunt=1;

%НАСТРОЙКА КОНТУРА СКОРОСТИ (симметричный оптимум)
Tw=8*T;                         % - speed regulator time constant
Kfb_w=1/wn;                     % - speed feedback coefficient
Kw=Kfb_i*(1/Cma)*Js/(4*Tconv*Kfb_w); % - speed regulator gain in pu
% Kw=(Js/(4*Tconv))              % - speed regulator gain

% Limitations
Il=10*In;                        % - Current limit value, A; 
Isc=Udc/(Ra*3/2);

kl=Il/Isc;                      % - Current limitation coefficient;
kli=Il*Kfb_i;                   % - current limitation value;

Klregs=1;                       % - speed regulator limitation;
Klregi=kli;                     % - current regulator limitation;
Klship=1;                       % - PWM input limitation;
