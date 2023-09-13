data = load('dataset.mat');
%% 
t = data.cigan.X.Data';
u = data.cigan.Y(:, 3).Data';
w = data.cigan.Y(:, 2).Data';
poz = data.cigan.Y(:, 1).Data';

plot(t, [200*u, w]); hold on;
xlabel('Timp');
ylabel('Amplitudine');

%% 

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

A = -1/T;
B = K/T; 
C = 1;
D = 0;
sysOne = ss(A,B,C,D);
ysimOne = lsim(sysOne, u_new-u_new(1), t, w(1));

plot(t, ysimOne); 
title('Identificare clasica - Neparametrica');
legend('Factor de umplere PWM', 'Viteza', 'Y63', 'Modelul');
xlabel('Timp');
ylabel('Amplitudine');

J = getRelativeError(w, ysimOne)
eMPN = getNormalizedMeanSquaredError(w, ysimOne)
%% viteza - pozitie neparametrica

i8 = 7434;
i9 = 7687;

K_poz = (poz(i9)-poz(i8))/(mean(w(i8:i9))*(t(i9)-t(i8)));

A1 = [-1/T 0; K_poz 0];
B1 = [K/T; 0];
C1 = [1 0 ; 0 1];
D1 = [0; 0];

sysTwo = ss(A1, B1, C1, D1);
ysimTwo = lsim(sysTwo , u, t, [w(1), poz(1)]);
ysimTwo = ysimTwo(:,2);

figure
plot(t, ysimTwo, t, poz), 
title("Suprapunerea intre integrata vitezei peste pozitie")
legend('Pozitia identificata', 'Pozitia propriu zisa')
xlabel('Timp');
ylabel('Amplitudine');

eMPN_pozitie = norm(poz - ysimTwo) / norm(poz - mean(poz));

%% Functia de transfer cu timp mort

Hs = tf(K, [T, 1], 'IODelay', Tau_m);

%% Preluare date din grafic

u = double(u);
t = double(t);
w = double(w);
poz = double(poz);

i1 = 1155;
i2 = 3058;
i3 = 4475;
i4 = 6208;
pas = 8;

t_id = u(i1:pas:i2);
u_id = u(i1:pas:i2);
w_id = w(i1:pas:i2);
p_id = poz(i1:pas:i2);

t_vd = u(i3:pas:i4);
u_vd = u(i3:pas:i4);
w_vd = w(i3:pas:i4);
p_vd = poz(i3:pas:i4);

Te = pas * (t(5) - t(4));
data_id = iddata(w_id, u_id, Te);
data_vd = iddata(w_vd, u_vd, Te);

%% Autocorelatie Intrare-Viteza cu ARMAX
sys_armax = armax(data_id, [1, 1, 1, 1]); % inclus gradul mort

% Scoatem functia de transfer
Hw_armax_viteza = tf(sys_armax.B, sys_armax.A, Te, 'variable', 'z^-1');
Hw_armax_viteza_continuu = d2c(Hw_armax_viteza, 'zoh');

K = 3172 / 12.97;
T = 1 / 12.97;

Hfinal = tf(K, [T, 1]);

A = -1/T;
B = K/T; 
C = 1;
D = 0;

sysOne = ss(A,B,C,D);
ysim = lsim(sysOne, u_new, t, w(1));

eMPN = norm(w-ysim)/norm(w-mean(w));
fprintf(eMPN);

plot(t, [w, ysim]);
legend('Viteza propriu zisa', 'Viteza identificata')
xlabel('Timp');
ylabel('Amplitudine');

% Plotare
figure;
subplot(1, 2, 1)
compare(data_vd, sys_armax)
title('Compare Viteza-Intrare cu ARMAX');
subplot(1, 2, 2)
resid(data_vd, sys_armax)
title('Resid Viteza-Intrare cu ARMAX');

%% Intercorelatie / XCorelatie Intrare-Viteza cu iv4

sys_iv4 = iv4(data_id, [1, 1, 1]);

% Scoatem functia de transfer
Hw_iv4_viteza = tf(sys_iv4.B, sys_iv4.A, Te, 'variable', 'z^-1');
Hw_iv4_viteza_continuu = d2c(Hw_iv4_viteza, 'zoh');

K = 5210 / 21.4;
T = 1 / 21.4;
HFinal2 = tf(K, [T, 1]);

A = -1/T;
B = K/T; 
C = 1;
D = 0;

sysTwo = ss(A,B,C,D);
ysim = lsim(sysTwo, u_new, t, w(1));

eMPN = norm(w-ysim)/norm(w-mean(w));
fprintf(eMPN);

plot(t, [w, ysim]);
legend('Viteza propriu zisa', 'Viteza identificata')
xlabel('Timp');
ylabel('Amplitudine');

% Plotare
figure;
subplot(1, 2, 1)
compare(data_vd, sys_iv4)
title('Compare intrare-viteza cu iv4');
subplot(1, 2, 2)
resid(data_vd, sys_iv4)
title('Resid intrare-viteza cu iv4');

%% Autocorelatie Viteza-Pozitie cu ARX

data_id_poz_vit = iddata(p_id, w_id, Te);
data_vd_poz_vit = iddata(p_vd, w_vd, Te);

sys_armax_poz_vit = arx(data_id_poz_vit, [1, 1, 0]);

% Scoatem functia de transfer
Hw_armax_pozitie_viteza = tf(sys_armax_poz_vit.B, sys_armax_poz_vit.A, Te, 'variable', 'z^-1');
Hw_armax_pozitie_viteza_continuu = tf(0.02841 / Te, [1, 0]);

fprintf(Hw_armax_pozitie_viteza_continuu)

A = 0;
B = 4.932;
C = 1;
D = 0;

sysThree = ss(A,B,C,D);
ysim = lsim(sysThree, w, t, poz(1));

plot(t, [poz, ysim]);
legend('Pozitia propriu zisa', 'Pozitia identificata')
xlabel('Timp');
ylabel('Amplitudine');

eMPN = norm(poz-ysim)/norm(poz-mean(poz));
fprintf(eMPN);

% Plotare
figure;
subplot(1, 2, 1)
compare(data_vd_poz_vit, sys_armax_poz_vit)
title('Compare Pozitie-Viteza cu ARMAX');
subplot(1, 2, 2)
resid(data_vd_poz_vit, sys_armax_poz_vit)
title('Resid Pozitie-Viteza cu ARMAX');

%% Intercorelatie / XCorelatie Viteza-Pozitie cu OE

data_id_poz_vit = iddata(p_id, w_id, Te);
data_vd_poz_vit = iddata(p_vd, w_vd, Te);
sys_oe_poz_vit = oe(data_id_poz_vit, [1, 1, 0]); % inclus gradul mort

% Scoatem functia de transfer
Hw_oe_pozitie_viteza = tf(sys_oe_poz_vit.B, sys_oe_poz_vit.F, Te, 'variable', 'z^-1');
Hw_armax_pozitie_viteza_continuu = tf(0.02841 / Te, [1, 0]);

A = 0;
B = 4.932;
C = 1;
D = 0;

sysThree = ss(A,B,C,D);
ysim = lsim(sysThree, w, t, poz(1));

plot(t, [poz, ysim]);
legend('Pozitia propriu zisa', 'Pozitia identificata')
xlabel('Timp');
ylabel('Amplitudine');

eMPN = norm(poz-ysim)/norm(poz-mean(poz));

% Plotare
figure;
subplot(1, 2, 1)
compare(data_vd_poz_vit, sys_oe_poz_vit)
title('Compare Pozitie-Viteza cu OE');
subplot(1, 2, 2)
resid(data_vd_poz_vit, sys_oe_poz_vit)
title('Resid Pozitie-Viteza cu OE');