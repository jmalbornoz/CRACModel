function [Pcrac f] = crac_model_v1(Q, Ts, Tr, P, phi, Pmax, fmax, e_coil, e_fan)
%%
%% crac_model_v1.m
%%
%% CRAC model
%% Inputs to the model:
%%   * Q: amount of heat to be removed, watss
%%   * Ts: supply air temperature, Celsius
%%   * Tr: return air temperature, Celsius
%%   * P: atmospheric pressure, pascals
%%   * phi: relative humidity, %
%%   * Pmax: maximum power consumed by the CRAC fans, watts
%%   * fmax: maximum airflow rate delivered by the CRAC fans, m^3/s
%%   * e_coil: cooling coil efficiency
%%   * e_fan: fan efficiency
%%
%%
%% Output:
%%   * Pcrac: power consumed by the CRAC fans
%%   * f: airflow rate, m^3/s
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%

%% Constants

% specific gas constant for dry air in Joule/(kg.K)
Rdry = 287.058;  

% specific gas constant for water vapour in Joule/(kg.K)
Rvap = 461.495;

% specific heat capacity of air between 0 and 40 degrees Celsius in Joule/(kg.K)
c_a = 1005;

% polymonial fit for specific heat capacity of water vapor
poly = [1.1255e-006 -4.2305e-004 1.8901];

%%

% water vapour saturation pressure in pascals
% temperatures in Celsius
Psat = 611.21*exp((18.678 - Ts/234.5).*Ts/(257.14 + Ts));

% water vapour partial pressure, pascals
Pvap = phi*Psat;

% dry air partial pressure, pascals
Pdry = P - Pvap;

% density of humid air, kg/m^3
rho = Pdry./(Rdry*(Ts + 273)) + Pvap./(Rvap*(Ts + 273));

% humidity ratio, kg H20 / kg dry air
H = 0.62198*Pvap./Pdry;

% specific heat capacity as a function of water vapour content, Joule/(kg.K)
c_ah = c_a + 1820*H; 

% calculates required airflow rate (m^3/s)
f = Q./(rho.*c_ah.*(Tr - Ts)*e_coil);

% calculates required fan power
%if f <= fmax
    Pcrac = Pmax*(f/fmax).^3;
    Pcrac = Pcrac/e_fan;
%else
%    warning('crac:maxflow','Required airflow rate exceeds maximum rated value')
%end

