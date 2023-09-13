data = load('dataset.mat');
%% 
t = data.cigan.X.Data';
u = data.cigan.Y(:, 3).Data';
w = data.cigan.Y(:, 2).Data';
pos = data.cigan.Y(:, 1).Data';

plot(t, [200*u, w]); hold on;
xlabel('Timp');
ylabel('Amplitudine');

%% intrare-viteza neparametric

i3 = 7434;
i4 = 7687;

ust = mean(u(i3:i4));
yst = mean(w(i3:i4));

K = yst / ust;
y63 = 0.63*yst;

plot(t, y63*ones(1, length(t)), 'r');
xlabel('Timp');
ylabel('Amplitudine');

i5 = 7189;
i6 = 7242;
T = t(i6) - t(i5);

Tm_i7 = 7151;
Tm_i8 = 7188; 

Tau_m = t(Tm_i8) - t(Tm_i7);

nidx = round(Tau_m / (t(12)-t(11)));
u_new = [u(1)*ones(nidx,1); u(1:length(u)-nidx)];

[A, B, C, D] = getStateSpaceForSysOrderOne(K, T);

sysUW = ss(A,B,C,D);
ysimUW = lsim(sysUW, u_new-u_new(1), t, w(1));

plot(t, ysimUW); 
title('Identificare clasica - Neparametrica');
legend('Factor de umplere PWM', 'Viteza', 'Y63', 'Modelul');
xlabel('Timp');
ylabel('Amplitudine');

J = getRelativeError(w, ysimUW)
eMPN = getNormalizedMeanSquaredError(w, ysimUW)
%% viteza - pozitie neparametrica

i8 = 7434;
i9 = 7687;

K_poz = (pos(i9)-pos(i8))/(mean(w(i8:i9))*(t(i9)-t(i8)));

A_pos = [-1/T 0; K_poz 0];
B_pos = [K/T; 0];
C_pos = [1 0 ; 0 1];
D_pos = [0; 0];

sysWPos = ss(A_pos, B_pos, C_pos, D_pos);
ysimWPos = lsim(sysWPos , u, t, [w(1), pos(1)]);
ysimWPos = ysimWPos(:,2);

figure
plot(t, ysimWPos, t, pos), 
title("Suprapunerea intre integrata vitezei peste pozitie")
legend('Pozitia identificata', 'Pozitia propriu zisa')
xlabel('Timp');
ylabel('Amplitudine');

eMPN_pozitie = getNormalizedMeanSquaredError(pos, ysimWPos)

%% Functia de transfer cu timp mort

Hs = tf(K, [T, 1], 'IODelay', Tau_m)

%% Preluare date din grafic si preprocesare

u = double(u);
t = double(t);
w = double(w);
pos = double(pos);

i1 = 1155;
i2 = 3058;
i3 = 4475;
i4 = 6208;
step = 8;

t_identification = u(i1:step:i2);
u_identification = u(i1:step:i2);
w_identification = w(i1:step:i2);
pos_identification = pos(i1:step:i2);

t_validation = u(i3:step:i4);
u_validation = u(i3:step:i4);
w_validation = w(i3:step:i4);
pos_validation = pos(i3:step:i4);

Te = step * (t(5) - t(4));
data_identification = iddata(w_identification, u_identification, Te);
data_validation = iddata(w_validation, u_validation, Te);

%% Autocorelatie Intrare-Viteza cu ARMAX
sys_armax = armax(data_identification, [1, 1, 1, 1]);

Hw_armax_w = tf(sys_armax.B, sys_armax.A, Te, 'variable', 'z^-1');
Hw_armax_w_continuous = d2c(Hw_armax_w, 'zoh');

% extracted from continuous Hw_armax_w_continuous
K = 3172 / 12.97;
T = 1 / 12.97;

Hfinal = tf(K, [T, 1]);

[A, B, C, D] = getStateSpaceForSysOrderOne(K, T);

sysArmaxUW = ss(A,B,C,D);
ysimUW = lsim(sysArmaxUW, u_new, t, w(1));

eMPN = getNormalizedMeanSquaredError(w, ysimUW)

plot(t, [w, ysimUW]);
legend('Viteza propriu zisa', 'Viteza identificata')
xlabel('Timp');
ylabel('Amplitudine');

figure;
subplot(1, 2, 1)
compare(data_validation, sys_armax)
title('Compare Viteza-Intrare cu ARMAX');
subplot(1, 2, 2)
resid(data_validation, sys_armax)
title('Resid Viteza-Intrare cu ARMAX');

%% Intercorelatie / XCorelatie Intrare-Viteza cu iv4

sys_iv4 = iv4(data_identification, [1, 1, 1]);

Hw_iv4_w = tf(sys_iv4.B, sys_iv4.A, Te, 'variable', 'z^-1');
Hw_iv4_w_continuous = d2c(Hw_iv4_w, 'zoh');

% extracted from continuous Hw_iv4_w_continuous
K = 5210 / 21.4;
T = 1 / 21.4;
HFinal2 = tf(K, [T, 1]);

[A, B, C, D] = getStateSpaceForSysOrderOne(K, T);

sysIv4UW = ss(A, B, C, D);
ysimUW = lsim(sysIv4UW, u_new, t, w(1));

eMPN = getNormalizedMeanSquaredError(w, ysimUW)

plot(t, [w, ysimUW]);
legend('Viteza propriu zisa', 'Viteza identificata')
xlabel('Timp');
ylabel('Amplitudine');

figure;
subplot(1, 2, 1)
compare(data_validation, sys_iv4)
title('Compare intrare-viteza cu iv4');
subplot(1, 2, 2)
resid(data_validation, sys_iv4)
title('Resid intrare-viteza cu iv4');

%% Autocorelatie Viteza-Pozitie cu ARX

data_id_pos_w = iddata(pos_identification, w_identification, Te);
data_vd_pos_w = iddata(pos_validation, w_validation, Te);

sys_armax_pos_w = arx(data_id_pos_w, [1, 1, 0]);

Hw_arx_pos_w = tf(sys_armax_pos_w.B, sys_armax_pos_w.A, Te, 'variable', 'z^-1');
Hw_arx_pos_w_continuous = tf(0.02841 / Te, [1, 0])

% extracted from continuous Hw_arx_pos_w_continuous
A = 0;
B = 4.932;
C = 1;
D = 0;

sysArxPosW = ss(A,B,C,D);
ysimUW = lsim(sysArxPosW, w, t, pos(1));

plot(t, [pos, ysimUW]);
legend('Pozitia propriu zisa', 'Pozitia identificata')
xlabel('Timp');
ylabel('Amplitudine');

eMPN = getNormalizedMeanSquaredError(pos, ysimUW)

figure;
subplot(1, 2, 1)
compare(data_vd_pos_w, sys_armax_pos_w)
title('Compare Pozitie-Viteza cu ARMAX');
subplot(1, 2, 2)
resid(data_vd_pos_w, sys_armax_pos_w)
title('Resid Pozitie-Viteza cu ARMAX');

%% Intercorelatie / XCorelatie Viteza-Pozitie cu OE

data_id_pos_w = iddata(pos_identification, w_identification, Te);
data_vd_pos_w = iddata(pos_validation, w_validation, Te);
sys_oe_pos_w = oe(data_id_pos_w, [1, 1, 0]); % inclus gradul mort

Hw_oe_pos_w = tf(sys_oe_pos_w.B, sys_oe_pos_w.F, Te, 'variable', 'z^-1');
Hw_oe_pos_w_continuous = tf(0.02841 / Te, [1, 0]);

% extracted from continuous Hw_oe_pos_w_continuous
A = 0;
B = 4.932;
C = 1;
D = 0;

sysOePosW = ss(A,B,C,D);
ysimUW = lsim(sysOePosW, w, t, pos(1));

plot(t, [pos, ysimUW]);
legend('Pozitia propriu zisa', 'Pozitia identificata')
xlabel('Timp');
ylabel('Amplitudine');

eMPN = getNormalizedMeanSquaredError(pos, ysimUW)

figure;
subplot(1, 2, 1)
compare(data_vd_pos_w, sys_oe_pos_w)
title('Compare Pozitie-Viteza cu OE');
subplot(1, 2, 2)
resid(data_vd_pos_w, sys_oe_pos_w)
title('Resid Pozitie-Viteza cu OE');