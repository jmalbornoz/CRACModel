%%
%% density and specific heat capacity of dry air
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%
clear all
close all

load -ascii 'air_properties.dat'

T = air_properties(:,1);
T = T(2:length(T));
rho = air_properties(:,2);
rho = rho(2:length(rho));
c = air_properties(:,3);
c = c(2:length(c));


Rdry = 287.058     % specific gas constant for dry air in Joule/(kg . K)
P = 101325;    % pressure in Pascals
Ta = 0;     % air temperature in Celsius

% formula for density of dry air
rho_exp = P./(Rdry*(T + 273));

figure(1)
plot(T,rho,'o',T,rho_exp,'r')
grid
xlabel('Temperature, Celsius')
ylabel('Density, kg/m^3')
legend('Data','Formula')
title('Density of dry air')

figure(2)
plot(T,c)
grid
xlabel('Temperature, Celsius')
ylabel('Specific heat capacity, kJ/(kg.K)')
title('Specific heat capacity of dry air')


%%
%% Humid air
%%

% range of data centre temperatures
T = linspace(5,40);

% density of dry air
rho_dry = P./(Rdry*(T + 273));

% specific gas constant for water vapour in Joule/(kg . K)
Rvap = 461.495;

% relative humidity 
phi = 0.4;

% water vapour saturation pressure in Pa, temperatures in Celsius
Psat = 611.21*exp((18.678 - T./234.5).*T./(257.14 + T));

% water vapour partial pressure
Pvap = phi*Psat;

% dry air partial pressure
Pdry = P - Pvap;

% density of humid air
rho_humid = Pdry./(Rdry*(T + 273)) + Pvap./(Rvap*(T + 273));

figure(3)
plot(T,rho_dry,T,rho_humid)
grid
xlabel('Temperature, Celsius')
ylabel('Density, kg/m^3')
legend('Dry air','Humid air')
title('Densities of dry and humid air')

% percent difference in dry and moist air densities
rho_diff = (rho_dry - rho_humid)*100./rho_dry;

figure(4)
plot(T,rho_diff)
grid
xlabel('Temperature, Celsius')
ylabel('Percentage difference')
title('Impact of humidity on air density')