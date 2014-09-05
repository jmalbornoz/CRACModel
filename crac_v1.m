%%
%% crac_v1.m
%%
%% CRAC model
%% Inputs to the model:
%%   * Q: amount of heat to be removed
%%   * Ts: air supply temperature
%%   * Tr: air return temperature
%%   * P: atmospheric pressure
%%   * phi: relative humidity
%%   * e: cooling coil efficiency
%%   * Pmax: maximum power consumed by the CRAC fans
%%   * fmax: maximum airflow delivered by the CRAC fans
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%
clear all
close all

%% Constants

% specific gas constant for dry air in Joule/(kg . K)
Rdry = 287.058;  

% specific gas constant for water vapour in Joule/(kg . K)
Rvap = 461.495;

% specific heat capacity of air between 0 and 40 degrees Celsius in J/(kg.K)
c_a = 1005;

%%
%% Case 1: dry air
%% 

% cooling coil efficiency
e = 1.0;

% maximum power consumed by CRAC fans, W
Pmax = 10e3;

% maximum airflow delivered by the CRAC fans, m^3/s
fmax = 6;

% atmospheric pressure, Pascals
P = 101325;

% amount of heat to be removed, Watts
Q = 100e3;

% air supply temperature, Celsius
Ts = 20;

% air return temperature, Celsius
Tr = 35;

Tair = 5:40;

% calculates air density at the supply temperature
rho =  P./(Rdry*(Tair + 273));

% calculates required airflow rate (m^3/s)
f_dry = Q./(rho*c_a*(Tr - Ts)*e);

% calculates required fan power
if f_dry <= fmax
    Pcrac_dry = Pmax*(f_dry/fmax).^3; 
else
    error('crac:maxflow','Required airflow rate exceeds maximum rated value')
end

figure(1)
plot(Tair,Pcrac_dry/1e3)
grid
xlabel('Air temperature, Celsius')
ylabel('CRAC power, kW')
title('CRAC power consumption as a function of temperature, dry air')

%%
%% Case 2: humid air
%%

% relative humidity 
phi = 0.55;

% water vapour saturation pressure in Pa, temperatures in Celsius
Psat1 = 611.21*exp((18.678 - Tair./234.5).*Tair./(257.14 + Tair));
Psat2 = 610.78*10.^((7.5*(Tair + 273) - 2048.625)./(Tair + 237.15));

figure(2)
plot(Tair,Psat1,Tair,Psat2)
grid
xlabel('Air temperature, Celsius')
ylabel('Water saturation pressure, Pascals')
title('Water vapour saturation pressure')

% water vapour partial pressure
Pvap = phi*Psat1;

% dry air partial pressure
Pdry = P - Pvap;

% specific humidity
H = 0.62198*Pvap/Pdry;

% density of humid air
rho_humid = Pdry./(Rdry*(Tair + 273)) + Pvap./(Rvap*(Tair + 273));

% calculates required airflow rate (m^3/s), constant specific heat value
f_humid1 = Q./(rho_humid*c_a*(Tr - Ts)*e);

% calculates required fan power
if f_humid1 <= fmax
    Pcrac_humid1 = Pmax*(f_humid1/fmax).^3; 
else
    error('crac:maxflow','Required airflow rate exceeds maximum rated value')
end

% the variation in specific heat capacity of humid air is now considered:

% specific heat capacity as a function of water vapour content
c_aw = c_a + 1820*H; 

% calculates required airflow rate (m^3/s), specific heat changes with
% temperature
f_humid2 = Q./(rho_humid*c_aw*(Tr - Ts)*e);

% calculates required fan power
if f_humid2 <= fmax
    Pcrac_humid2 = Pmax*(f_humid2/fmax).^3; 
else
    error('crac:maxflow','Required airflow rate exceeds maximum rated value')
end

figure(3)
plot(Tair,Pcrac_dry/1e3,Tair,Pcrac_humid1/1e3,Tair,Pcrac_humid2/1e3)
grid
xlabel('Temperature, Celsius')
ylabel('CRAC power, kW')
legend('Dry air','Humid air (c_a)','Humid air (c_aw)')
title('CRAC power consumption as a function of temperature')

% the influence of changing specific heat capacity of air with humidity is
% negligible

%%
%% Now we explore constant temperature and changing humidity
%%

% air supply temperature
Ts = 20;

% relative humidity within ASHRAE's limits
phi = linspace(0.4,0.65);

% dry air

% calculates air density at the supply temperature
rho =  P./(Rdry*(Ts + 273));

% calculates required airflow rate (m^3/s)
f_dry = Q./(rho*c_a*(Tr - Ts)*e);

% calculates required fan power
if f_dry <= fmax
    Pcrac_dry = Pmax*(f_dry/fmax).^3;
    Pcrac_dry = Pcrac_dry*ones(1,length(phi));
else
    error('crac:maxflow','Required airflow rate exceeds maximum rated value')
end

% humid air

% water vapour saturation pressure in Pa, temperatures in Celsius
Psat = 611.21*exp((18.678 - Ts./234.5).*Ts./(257.14 + Ts));

% water vapour partial pressure
Pvap = phi*Psat;

% dry air partial pressure
Pdry = P - Pvap;

% specific humidity
H = 0.62198*Pvap./Pdry;

% specific heat capacity as a function of water vapour content
c_aw = c_a + 1820*H; 

% density of humid air
rho_humid = Pdry./(Rdry*(Ts + 273)) + Pvap./(Rvap*(Ts + 273));

% calculates required airflow rate (m^3/s)
f_humid = Q./(rho_humid.*c_aw*(Tr - Ts)*e);

% calculates required fan power
if f_humid <= fmax
    Pcrac_humid = Pmax*(f_humid/fmax).^3; 
else
    error('crac:maxflow','Required airflow rate exceeds maximum rated value')
end

figure(4)
plot(phi,Pcrac_dry/1e3,phi,Pcrac_humid/1e3)
grid
xlabel('Relative humidity')
ylabel('CRAC power, kW')
legend('Dry air','Humid air')
title('CRAC power consumption as a function of relative humidity')
