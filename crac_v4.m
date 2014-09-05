%%
%% crac_v3.m
%%
%% This file uses the CRAC model defined in crac_model_v2.m
%% Inputs:
%%   * Q: amount of heat to be removed in watts
%%   * Ts: supply air temperature, Celsius
%%   * Tr: return air temperature, Celsius
%%   * Tcoil: cooling coil temperature, Celsius
%%   * P: atmospheric pressure in pascals
%%   * phi: relative humidity, %
%%   * Pmax: maximum power consumed by the CRAC fans in watts
%%   * fmax: maximum airflow delivered by the CRAC fans in m^3/s
%%   * e: cooling coil efficiency
%%   * mode: ASHRAE's humidity range to be adopted
%%           mode = 1: 2004 recommended range
%%           mode = 2: 2004 allowable range
%%           mode = 3: 2008 recommended range
%%           mode = 4: 2008 allowable range
%%
%%   * change in specific enthalpy of air is used to calculate requiredd
%%     airflow
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%
clear all
close all

% amount of heat to be removed, Watts
Q = 30e3;

% air supply temperature, Celsius
Ts = 18;

% air return temperature, Celsius
Tr = 30;

% cooling coil temperature = entering water temperature, Celsius
Tcoil = 7.22;

% atmospheric pressure, Pascals
P = 101325;

% relative humidity
%phi = 0.5;

% maximum rated power for CRAC fans, W
Pmax = 2.2e3;

% maximum airflow delivered by the CRAC fans, m^3/s
fmax = 2.38;

% cooling coil efficiency
e_coil = 1.0;

% fan efficency
e_fan = 1.0;

%%
%% Variation of fan power consumption with changing humidity
%%

%relative humidity
phi = linspace(0.4,0.55);
%phi = 0.323;

P_max = Pmax*ones(1,length(phi));
f_max = fmax*ones(1,length(phi));

% fan power draw, W
Pcrac = zeros(1,length(phi));
f = zeros(1,length(phi));
phi_s = zeros(1,length(phi));
Qload = zeros(1,length(phi));
for i = 1: length(phi)
    [Pcrac(i) f(i) phi_s(i) Qload(i)] = crac_model_v3(Q, Ts, Tr, Tcoil, P, phi(i), Pmax, fmax, e_coil, e_fan, 1);
end
    
figure(1)
%plot(phi,Pcrac/1e3,phi,P_max/1e3)
plot(phi,Pcrac/1e3)
grid
xlabel('Relative humidity')
ylabel('Power draw (kW)')
%title('Fan power draw as a function of relative humidity')
%legend('Fan power draw','Maximum fan rated power')
axis([0.4 0.55 1.356 1.37])

figure(2)
%plot(phi,f,phi,f_max)
plot(phi,Qload/1e3)
grid
xlabel('Relative humidity')
ylabel('Total heat output (kW)')
%title('Required airflow as a function of relative humidity')
axis([0.4 0.55 42 57])
%legend('Airflow','Maximum rated airflow')

figure(3)
plot(phi,phi_s*100)
grid
xlabel('Air supply temperature')
ylabel('Relative humidity (%)')
%title('Required airflow as a function of supply temperature')

%%
%% Variation of fan power consumption with changing temperature
%%

%relative humidity
phi = 0.55;

% air supply and return temperatures, Celsius
Ts = linspace(20,25);
deltaT = 11;
Tr = linspace(20 + deltaT,25 + deltaT);
%Tr = 36;

% fan power draw, W
Pcrac = zeros(length(Ts));
f = zeros(length(Tr));
phi_s = zeros(length(Ts));
Qload = zeros(length(Ts));
for i = 1: length(Ts)
    [Pcrac(i) f(i) phi_s(i) Qload(i)] = crac_model_v3(Q, Ts(i), Tr(i), Tcoil, P, phi, Pmax, fmax, e_coil, e_fan, 1);
    %[Pcrac(i) f(i) phi_s(i) Qload(i)] = crac_model_v3(Q, Ts(i), Tr, Tcoil, P, phi, Pmax, fmax, e_coil, e_fan, 1);
end

figure(4)
plot(Ts,Pcrac/1e3)
grid
xlabel('Air supply temperature (Celsius)')
ylabel('Power draw (kW)')
axis([20 25 1.6 2])
%title('CRAC power draw as a function of supply temperature')

figure(5)
plot(Ts,phi_s*100)
grid
xlabel('Air supply temperature (Celsius)')
ylabel('Relative humidity (%)')
axis([20 25 70 80])
%title('Required airflow as a function of supply temperature')

figure(6)
plot(Ts,Qload/1e3)
grid
xlabel('Air supply temperature (Celsius)')
ylabel('Total heat output (kW)')
%title('Required airflow as a function of relative humidity')
%axis([0.4 0.55 42 57])
%legend('Airflow','Maximum rated airflow')

%%
%% Simultaneous variation of supply humidity and temperature
%%

% air supply and return temperatures, Celsius
Ts = linspace(20,25);
deltaT = 11;
Tr = linspace(20 + deltaT,25 + deltaT);

%relative humidity
phi = linspace(0.4,0.55);

Pcrac = zeros(length(Ts),length(phi));
f = zeros(length(Ts),length(phi));
phi_s = zeros(length(Ts),length(phi));
Qload = zeros(length(Ts),length(phi));

for i = 1:length(Ts)
    for j = 1:length(phi)
            [Pcrac(i,j) f(i,j) phi_s(i,j) Qload(i,j)] = crac_model_v3(Q, Ts(i), Tr(i), Tcoil, P, phi(j), Pmax, fmax, e_coil, e_fan, 1);
    end
end

figure(7)
colormap(autumn)
mesh(phi*100,Ts,Pcrac/1e3)
xlabel('Return air relative humidity (%)')
ylabel('Supply air temperature (Celsius)')
zlabel('Fan power draw (kW)')
axis([40 55 20 25 1.78 1.90])

figure(8)
colormap(autumn)
mesh(phi*100,Ts,Qload/1e3)
axis([40 55 20 25  45 70])
xlabel('Return air relative humidity (%)')
ylabel('Supply air temperature (Celsius)')
zlabel('Total heat output (kW)')
