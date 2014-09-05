%%
%% Dry cooler load-temperature model
%%
%% * Affinity law has been considered to build this model.
%% * Consumption data from BAC drycooler model S8013-S213B (123 kW)
%%   3 fans, each driven by a 6-pole motor rated at 1.5 kW 
%%
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% June 2011
%%
clear all
close all

Qmax = 123e3;    % rated duty, W
EWT = 40;        % entering water temperature in cooling coils, celsius
P_fans = 4.5e3;  % rated fan consumption, W

deltaT_min = 12;              % minimum temperature difference in dry cooler, celsius
Tmax = EWT - deltaT_min;      % maximum temperature for dry cooler operation

deltaT = 1;

L = (0:0.05:1);
Q = L*Qmax;
T = -10:deltaT:8;

% calculates consumed power
pwr = zeros(length(T),length(Q));
for i = 1:length(T)
    for j = 1:length(Q)
        pwr(i,j) = P_fans*((Q(j)/(EWT - T(i)))/(Qmax/deltaT_min))^3;
    end
end

figure(1)
colormap(autumn)
mesh(L,T,pwr/1e3)
%surf(Q/1e3,T,pwr/1e3)
%title(name)
xlabel('Heat load (kW)')
ylabel('Dry bulb temperature (celsius)')
zlabel('Consumed power (kW)')
%legend('Data','Fit')
%hold off

pwr = pwr/Qmax;