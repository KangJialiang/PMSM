clc
clear
close all;
%% 100kW PMSM

syms i_d i_q t_e omega_r %设定变量

R_s = 8.3e-3; %定子相电阻
psi_f = 0.071; %定子匝链相绕组永磁磁链峰值
p = 4; %电机极对数
J = 0.1; %转动惯量
L_sigma = 30e-6; %定子绕组相漏电感
L_md = 174e-6; %d轴励磁电感
L_mq = 293e-6; %q轴励磁电感
n = 4700; %额定转速
T_emax = 256; %峰值转矩

L_d = L_sigma + L_md; %直轴同步电感
L_q = L_sigma + L_mq; %交轴同步电感

n_max = 12000; %最高转速

name = "100kW PMSM";

%% 公式

constT(i_d, i_q, t_e) = i_q == 2 * t_e / (3 * p * (psi_f + (L_d - L_q) * i_d)); %恒转矩电流表达式(5-61)
MTPA(i_d, i_q) = i_q == sqrt(i_d * (i_d + psi_f / (L_d - L_q))); %第二象限的MTPA曲线(5-64)
MTPA_sub(i_d, i_q) = i_q == -sqrt(i_d * (i_d + psi_f / (L_d - L_q))); %第三象限的MTPA曲线(5-64)
Te(i_d, i_q) = 3/2 * p * (psi_f * i_q + (L_d - L_q) * i_d * i_q); %转矩方程(5-53)
Pe(i_d, i_q, omega_r) = Te * omega_r / p; %输出功率(5-86)
Is(i_d, i_q) = sqrt(i_d^2 + i_q^2); %电流
Us(i_d, i_q, omega_r) = sqrt((omega_r * L_d * i_d + omega_r * psi_f)^2 + (omega_r * L_q * i_q)^2); %电压（忽略铜损）
MTPV(i_d, i_q) = i_q == sqrt(L_d * (L_d * i_d + psi_f) * ((L_d - L_q) * i_d + psi_f) / (L_d - L_q) / L_q^2); %第二象限的MTPV曲线(5-89)

%% 1）峰值转矩对应峰值相电流

[i_dmax, i_qmax] = vpasolve([constT(i_d, i_q, T_emax), MTPA], [i_d, i_q], [-inf, 0; 0, inf]); %数值求解两曲线交点
i_dmax = double(i_dmax);
i_qmax = double(i_qmax);
fprintf('峰值d轴电流：%f\t峰值q轴电流：%f\n', i_dmax, i_qmax);

Ismax = sqrt(i_dmax^2 + i_qmax^2); %最大定子电流
curLim(i_d, i_q) = i_d^2 + i_q^2 == Ismax^2; %电流极限圆(5-67)

%% 2）恒转矩曲线（渐近线）、等转矩线

i_dasymptote = psi_f / (L_q - L_d); %恒转矩曲线的竖直渐近线
i_qasymptote = 0; %恒转矩曲线的水平渐近线

figure(1);
fcontour(Te, [-Ismax, 2 * i_dasymptote, -Ismax, Ismax], 'LevelList', linspace(10, T_emax, 10)); %恒转矩曲线
c = colorbar;
c.Label.String = 'T_e/N\cdotm';
hold on
line([i_dasymptote, i_dasymptote], [-Ismax, Ismax], 'linestyle', '--'); %竖直渐近线
line([-Ismax, 2 * i_dasymptote], [i_qasymptote, i_qasymptote], 'linestyle', '--'); %水平渐近线
hold off
xlabel('$i_d$/A', 'Interpreter', 'latex')
ylabel('$i_q$/A', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + '恒转矩曲线（渐近线）')

figure(2);
fcontour(Te, [-Ismax, 0, -Ismax, Ismax], 'LevelList', linspace(-T_emax, T_emax, 20), 'LineStyle', '--'); %等转矩线
c = colorbar;
c.Label.String = 'T_e/N\cdotm';
hold on
fimplicit(MTPA, 'Color', 'b');
fimplicit(MTPA_sub, 'Color', 'b');
hold off
xlabel('$i_d$/A', 'Interpreter', 'latex')
ylabel('$i_q$/A', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + '等转矩线')

%% 3）电压极限椭圆（簇）

omega_rb = n * 2 * pi / 60 * p; %基速
omega_rmax = n_max * 2 * pi / 60 * p; %最高转速
Usm = omega_rb * sqrt((psi_f + L_d * i_dmax)^2 + (L_q * i_qmax)^2);

voltLim(i_d, i_q, omega_r) = (L_d * i_d + psi_f)^2 + (L_q * i_q)^2 == (Usm / omega_r)^2; %电压极限椭圆簇(5-74)

figure(3)
hold on

for omega = linspace(omega_rb, omega_rmax, 10)
    fimplicit(voltLim(i_d, i_q, omega), [-Ismax * 2.5, Ismax * 1.2, -Ismax, Ismax], 'LineStyle', '--'); %电压极限椭圆
end

fimplicit(curLim(i_d, i_q)); %电流极限圆
hold off
xlabel('$i_d$/A', 'Interpreter', 'latex')
ylabel('$i_q$/A', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + '电压极限椭圆和电流极限圆')

%% 4、5）最大转矩电流比（MTPA）曲线和最大转矩电压比（MTPV）曲线

figure(4)
hold on

for omega = linspace(omega_rb, omega_rmax, 10)
    fimplicit(voltLim(i_d, i_q, omega), [-Ismax * 2.5, Ismax * 1.2, -Ismax, Ismax], 'LineStyle', '--'); %电压极限椭圆
end

fimplicit(curLim(i_d, i_q)); %电流极限圆
fimplicit(MTPA, 'LineWidth', 1, 'MeshDensity', 300); %MTPA曲线
fimplicit(MTPV, 'LineWidth', 1, 'MeshDensity', 300); %MTPV曲线
hold off
xlabel('$i_d$/A', 'Interpreter', 'latex')
ylabel('$i_q$/A', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + ' MTPA曲线和MTPV曲线')

%% 6） 机械特性（外特性）

if psi_f / L_d > Ismax %图5-18a
    omega_r2 = omega_rmax;
else %图5-18b
    [i_d2, i_q2] = vpasolve([MTPV, curLim], [i_d, i_q], [-inf, 0; 0, inf]); %A2点
    i_d2 = double(i_d2);
    i_q2 = double(i_q2);
    omega_r2 = vpasolve(voltLim(i_d2, i_q2, omega_r), omega_r, [0, inf]);
    omega_r2 = double(omega_r2);
end

omegaList = [linspace(0, omega_rb, 10), linspace(omega_rb, omega_r2, 20), linspace(omega_r2, omega_rmax, 10)];
nList = omegaList / 2 / pi * 60 / p;
idList = zeros(1, length(omegaList));
iqList = zeros(1, length(omegaList));

for i = 1:length(omegaList)

    if omegaList(i) <= omega_rb %恒转矩区
        idList(i) = i_dmax;
        iqList(i) = i_qmax;
    elseif omegaList(i) <= omega_r2 %弱磁I区
        [idtmp, iqtmp] = vpasolve([voltLim(i_d, i_q, omegaList(i)), curLim(i_d, i_q)], ...
            [i_d, i_q], [-inf, 0; 0, inf]); %电压和电流同时达极限
        idList(i) = double(idtmp);
        iqList(i) = double(iqtmp);
    elseif omegaList(i) <= omega_rmax %弱磁II区
        [idtmp, iqtmp] = vpasolve([voltLim(i_d, i_q, omegaList(i)), MTPV], ...
            [i_d, i_q], [-inf, 0; 0, inf]); %MTPV线与电压极限椭圆交点
        idList(i) = double(idtmp);
        iqList(i) = double(iqtmp);
    end

end

figure(6)
hold on
plot(nList, Te(idList, iqList), 'DisplayName', 't_e/N\cdot m', 'LineWidth', 1)
plot(nList, Pe(idList, iqList, omegaList) / 1e3, 'DisplayName', 'P_e/kW', 'LineWidth', 1)
plot(nList, Us(idList, iqList, omegaList), 'DisplayName', 'u_s/V', 'LineWidth', 1)
plot(nList, Is(idList, iqList), 'DisplayName', 'i_s/A', 'LineWidth', 1)
line([omega_rb / 2 / pi * 60 / p, omega_rb / 2 / pi * 60 / p], [0, 500], 'linestyle', '--', 'Color', 'b', 'DisplayName', 'n_{rb}'); %n_rb
line([omega_r2 / 2 / pi * 60 / p, omega_r2 / 2 / pi * 60 / p], [0, 500], 'linestyle', '--', 'Color', 'g', 'DisplayName', 'n_{r2}'); %n_r2
legend('Location', 'best')
hold off
xlabel('$n/\rm{rpm}$', 'Interpreter', 'latex')
ylabel('$t_e,P_e,u_s,i_s$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
title(name + '恒转矩与恒功率运行（外特性）曲线')

%% 7）额定转速、最大转矩工作点的空间矢量图

figure(7)

hold on

scale = 1000; x = 0; y = 0; u = psi_f * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '{\psi}_f');
x = x + u; y = y + v;
u = 0; v = L_q * i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'L_qi_q');
x = x + u; y = y + v;
u = L_d * i_dmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'L_di_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\psi_s');

scale = 0.2; x = 0; y = 0; u = 0; v = i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_q');
x = x + u; y = y + v;
u = i_dmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_s');

scale = 1; x = 0; y = 0; u = 0; v = omega_rb * psi_f * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'e_0');
x = x + u; y = y + v;
u = R_s * i_dmax * scale; v = R_s * i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'R_si_s');
x = x + u; y = y + v;
u = -omega_rb * L_q * i_qmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\omega_sL_qi_q');
x = x + u; y = y + v;
u = 0; v = omega_rb * L_d * i_dmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\omega_sL_di_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'u_s');

legend('Location', 'bestoutside')
hold off
newcolors = {'#FFB6C1', '#FFC0CB', '#DC143C', '#FFF0F5', '#DB7093', '#FF69B4', ...
            '#FF1493', '#C71585', '#DA70D6', '#D8BFD8', '#DDA0DD', '#EE82EE'};
colororder(newcolors)
xlabel('$d$', 'Interpreter', 'latex')
ylabel('$q$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + '额定转速、最大转矩工作点的空间矢量图')

%% 8）额定转速、最大制动转矩回馈制动的空间矢量图

figure(8)

hold on

scale = 1000; x = 0; y = 0; u = psi_f * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '{\psi}_f');
x = x + u; y = y + v;
u = 0; v = L_q * -i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'L_qi_q');
x = x + u; y = y + v;
u = L_d * i_dmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'L_di_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\psi_s');

scale = 0.2; x = 0; y = 0; u = 0; v = -i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_q');
x = x + u; y = y + v;
u = i_dmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'i_s');

scale = 1; x = 0; y = 0; u = 0; v = omega_rb * psi_f * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'e_0');
x = x + u; y = y + v;
u = R_s * i_dmax * scale; v = R_s * -i_qmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'R_si_s');
x = x + u; y = y + v;
u = -omega_rb * L_q * -i_qmax * scale; v = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\omega_sL_qi_q');
x = x + u; y = y + v;
u = 0; v = omega_rb * L_d * i_dmax * scale;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', '\omega_sL_di_d');
u = x + u; v = y + v;
x = 0; y = 0;
quiver(x, y, u, v, "off", 'LineWidth', 2, 'DisplayName', 'u_s');

legend('Location', 'bestoutside')
hold off
newcolors = {'#FFB6C1', '#FFC0CB', '#DC143C', '#FFF0F5', '#DB7093', '#FF69B4', ...
            '#FF1493', '#C71585', '#DA70D6', '#D8BFD8', '#DDA0DD', '#EE82EE'};
colororder(newcolors)
xlabel('$d$', 'Interpreter', 'latex')
ylabel('$q$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + '额定转速、最大制动转矩回馈制动的空间矢量图')

%% 9）MTPA对应的最佳电流分配曲线

IsList = linspace(0, Ismax, 20);
idList = zeros(1, length(IsList));
iqList = zeros(1, length(IsList));

for i = 1:length(IsList)
    [idtmp, iqtmp] = vpasolve([MTPA, i_d^2 + i_q^2 == IsList(i)^2], [i_d, i_q], [-inf, 0; 0, inf]);
    idList(i) = double(idtmp);
    iqList(i) = double(iqtmp);
end

figure(9)
hold on
plot(IsList, idList, 'DisplayName', 'i_d', 'LineWidth', 1)
plot(IsList, iqList, 'DisplayName', 'i_q', 'LineWidth', 1)
legend('Location', 'northeast')
hold off
xlabel('$i_s/\rm{A}$', 'Interpreter', 'latex')
ylabel('$i_d,i_q/\rm{A}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis equal;
title(name + ' MTPA对应的最佳电流分配曲线')

%% 10）不同电流幅值下的"矩角特性曲线"

syms i_s beta

TeExci(i_s, beta) = 3/2 * p * psi_f * i_s * sin(beta); %励磁转矩
TeRelu(i_s, beta) = 3/2 * p * 1/2 * (L_d - L_q) * i_s^2 * sin(2 * beta); %磁阻转矩

figure(10)
hold on

for i = IsList
    fplot(TeExci(i, beta / 180 * pi) + TeRelu(i, beta / 180 * pi), [0, 180]); %角度单位为°
end

betaList = atan2(iqList, idList); %单位为rad
betaList(betaList == 0) = pi / 2; %corner situation
plot(betaList / pi * 180, TeExci(IsList, betaList) + TeRelu(IsList, betaList), 'LineStyle', '--'); %标出MTPA曲线
hold off
xticks(0:10:180)
xlabel('$\beta/\rm{^\circ}$', 'Interpreter', 'latex')
ylabel('$t_e/\rm{N\cdot m}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'XLim', [0, 180], 'YLim', [0, inf])
title(name + '不同电流幅值下的"矩角特性曲线"')

%% 11）磁阻转矩/电磁转矩曲线

figure(11)

subplot(2, 1, 1);
hold on

for i = IsList
    fplot(TeExci(i, beta / 180 * pi), [0, 180]);
end

hold off
xticks(0:10:180)
xlabel('$\beta/\rm{^\circ}$', 'Interpreter', 'latex')
ylabel('$t_e/\rm{N\cdot m}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'XLim', [0, 180], 'YLim', [0, inf])
title(name + '电磁转矩曲线')

subplot(2, 1, 2);
hold on

for i = IsList
    fplot(TeRelu(i, beta / 180 * pi), [0, 180]);
end

hold off
xticks(0:10:180)
xlabel('$\beta/\rm{^\circ}$', 'Interpreter', 'latex')
ylabel('$t_e/\rm{N\cdot m}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'XLim', [0, 180])
title(name + '磁阻转矩曲线')
