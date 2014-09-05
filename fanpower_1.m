%% fanpower_1.m
%%
%% Fan curves
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% December 2010
%%
clear all
close all

% loads file with system & fan curves
%
% flow_s, pressure_s: system curve 
% flow, pressure: fan curve
% flow_e, e: efficiency curve
%
load fan

% performs fits for system curve
P1 = polyfit(flow_s,pressure_s,2);

% system constant
Cs = 1/sqrt(P1(1));

% system curve
system = flow_s.^2/Cs^2;

% verifying system curve fit
figure(1)
plot(flow_s, pressure_s,'o')
hold on
plot(flow_s, system, 'r')
grid
xlabel('Air flow, m^3/s')
ylabel('Static pressure, Pascals')

% performs fit for fan curve
P2 = polyfit(flow,pressure,2);
fan = P2(1)*flow.^2 + P2(2)*flow + P2(3); 

% verifying fit for fan curve
plot(flow,pressure,'o')
plot(flow,fan,'g')
%legend('System data','System curve', ,'Fan curve')

% fit for efficiency curve
% a polynomial fit is no good
% we use a Fourier series approach instead
% the data is interpolated
flow = 0:0.1:23.2;
e_int = interp1(flow_e,e,flow);

% computes Fourier coefficients by numerical integration 
nterms = 50;
coeffs = [];

a0 = 0;
for i = 2:length(flow)
    sum = (flow(i) - flow(i - 1))*e_int(i);
    a0 = a0 + sum; 
end
a0 = a0/23.2;

for n = 1:nterms
    an = 0;
    mult = cos(2*n*pi*flow/23.2);
    g = e_int.*mult;
    for i = 2:length(flow)
        sum = (flow(i) - flow(i - 1))*g(i);
        an = an + sum; 
    end
    coeffs = [coeffs an];
end
coeffs = coeffs/23.2;

% builds Fourier series approximation
eff = zeros(1,length(flow));
a0 = a0*ones(1,length(flow));
for n = 1:nterms
    eff = eff + coeffs(n)*cos(2*n*flow*pi/23.2);
end
eff = a0 + 2*eff;

figure(2)
plot(flow_e,e,'o')
hold on
plot(flow,eff,'r')
grid
legend('Data','Fit')
%axis([0 24 0 1])

% re-defines flow, computes fan & system pressure
%flow = 0:1:24;
fan = P2(1)*flow.^2 + P2(2)*flow + P2(3);
system = flow.^2/Cs^2;

% computes power
power = flow.*fan;
power = power*4.7195e-4*249; 

figure(3)
plot(flow,system,flow,fan,'g',flow,power,'r')
grid
legend('System','Fan','Power')


%set(gca,'XTick',0:1:25)
%set(gca,'YTick',0:1:25)

%factor = 3.9;
%fit2 = factor*P1(1)*flow.^2 + factor*P1(2)*flow + P1(3);

%flow_OP = 14.4

%pressure_OP = factor*P1(1)*flow_OP.^2 + factor*P1(2)*flow_OP + P1(3)

%power = flow_OP*4.719e-4*pressure_OP*249

%flow = flow*4.7195e-4;      % flow in cubic meters per second
%pressure = pressure*249;    % pressure in Pascals

% test
%x = 0:0.01:pi;
%lim = 0;
%prev = rand(1)*0.3;
%x = [prev];
%while(prev < pi)
%    prev = prev + rand(1)*0.3;
%    x = [x prev];
%end

%y = sin(x);
%y = x.*ones(1,length(x));

% integral
%sum = 0;
%for i = 2:length(x)
%    sum = sum + (x(i) - x(i-1))*y(i);
%end
%figure(2)
%stem(flow_e,e)
%hold on
%plot(flow,e_int,'r')
%grid
