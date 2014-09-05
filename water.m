%%
%% specific heat capacity of water vapour
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%
clear all
close all

load -ascii 'water_vapour.dat'

% temperatures in Celsius
T = water_vapour(:,1) - 273;

% specific heat in joule/(kg-K)
c_wv = water_vapour(:,2)*1000;

% polynomial fit
poly = polyfit(T,c_wv,2);
fit = poly(1)*T.^2 + poly(2)*T + poly(3);

figure(1)
plot(T,c_wv,'*',T,fit)
grid
xlabel('Temperature, Celsius')
ylabel('Specific heat capacity, joule/(kg-K)')
title('Specific heat capacity of water vapour')

c_vap = @(z)(poly(1)*z.^2 + poly(2)*z + poly(3));

c_vap(0)
c_vap(35)



