function [Pcrac f phi_s Qload] = crac_model_v1(Q, Ts, Tr, Tcoil, P, phi, Pmax, fmax, e_coil, e_fan, mode)
%%
%% crac_model_v3.m
%%
%% CRAC model
%% Input:
%%   * Q: amount of heat to be removed in watts
%%   * Ts: supply air temperature in Celsius
%%   * Tr: return air temperature in Celsius
%%   * Tcoil: cooling coil temperature in Celsius
%%   * P: atmospheric pressure in pascals
%%   * phi: return air relative humidity, %
%%   * Pmax: maximum power consumed by the CRAC fans in watts
%%   * fmax: maximum airflow delivered by the CRAC fans in m^3/s
%%   * e_coil: cooling coil efficiency
%%   * e_fan: fan efficiency
%%   * mode: ASHRAE's humidity range to be adopted
%%           mode = 1: 2004 recommended range
%%           mode = 2: 2004 allowable range
%%           mode = 3: 2008 recommended range
%%           mode = 4: 2008 allowable range
%%
%% Output:
%%   * Pcrac: power consumed by the CRAC fans in watts
%%   * f: airflow rate in m^3/s
%%   * phi_s: supply relative humidity, %
%%   * Qload: heat load (sensible + latent) to next element in cooling chain, watts
%%
%% 21/11/2011: Dewpoint calculation is included 
%% 23/11/2011: Distinction between sensible and laten heat is included
%% 24/11/2011: c_a as a function of humidity for airflow calculation
%%
%% Jose Albornoz
%% Fujitsu Laboratories of Europe
%% November 2011
%%

%% Constants
%%

% specific gas constant for dry air in joule/(kg.K)
Rdry = 287.058;  

% specific gas constant for water vapour in joule/(kg.K)
Rvap = 461.495;

% specific heat capacity of air between 0 and 40 degrees Celsius in joule/(kg.C)
c_a = 1005;

% specific evaporation heat for water at 0 celsius in joule/kg
c_evap = 2501e3;

% molar mass of water, kg/mol 
M_H2O = 0.01801534;

% molar mass of dry air, kg/mol
M_dry = 0.0289644;

% specific heat capacity of water vapour, joule/(kg.C)
c_vap = 1860; 

%%
%%

% dewpoint calculation, Celsius
a = 17.271;
b = 237.7;
gamma = a*Tr/(b + Tr) + log(phi);
Tdew = b*gamma/(a - gamma);

% water vapour saturation pressure at return temperature
% pressures in Pascals, temperatures in Celsius
Psat = 611.21*exp((18.678 - Tr/234.5).*Tr/(257.14 + Tr));

% water vapour partial pressure at return temperature, pascals
Pvap = phi*Psat;

% dry air partial pressure, pascals
Pdry = P - Pvap;

% volume mixing ratio
X = Pvap/P;

% specific humidity in return air, kg H20/kg humid air 
H = X*M_H2O./(X*M_H2O + (1 - X)*M_dry);

% specific enthalpy of air at return temperature, joule/kg
hr = H*c_evap + H*Tr*c_vap + Tr*c_a;

%%

if Tcoil > Tdew     %% dry coil operation
    
    % specific enthalpy at supply temperature, joule/kg
    hs = H*c_evap + H*Ts*c_vap + Ts*c_a;
    
    % net enthalpy change, joule/kg
    delta_h = hr - hs;
    
    % water vapour saturation pressure after cooling, joule/kg
    Psat = 611.21*exp((18.678 - Ts/234.5).*Ts/(257.14 + Ts));
    
    % relative humidity after cooling
    phi_s = Pvap/Psat;   
    
    % air density after cooling, kg/m^3
    rho = Pdry./(Rdry*(Ts + 273)) + Pvap./(Rvap*(Ts + 273));

    % specific heat capacity of air as a function of water vapour content, Joule/(kg.C)
    c_ah = c_a + 1132*Pvap./Pdry;
    
    % required airflow for removal of sensible IT heat load, m^3/s
    f = Q./(rho.*c_ah*(Tr - Ts)*e_coil)
    
    %% total heat (sensible + latent), watts
    Qload = f*rho*delta_h;

    % calculates required fan power, watts
    Pcrac = Pmax*(f/fmax).^3;
    if Pcrac > Pmax
        error('Required CRAH power exceeds maximum rated value')
    end 
    Pcrac = Pcrac/e_fan;  
    Pdis = Pcrac*(1 - e_fan);
    
    % adds heat load + fan power draw
    Qload = Qload + Pdis;
 
    % message generation
    %switch mode
    %    case 1
    %        if phi_s < 0.4 || phi_s > 0.55
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2004 recommended range',phi_s*100,char(37))        
    %        end
    %    case 2
    %        if phi_s < 0.2 || phi_s > 0.8
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2004 allowable range',phi_s*100,char(37))        
    %        end         
    %    case 3
    %        % supply air dewpoint calculation, Celsius
    %        gamma = a*Ts/(b + Ts) + log(phi_s);
    %        Tdew_s = b*gamma/(a - gamma);
    %        if Tdew_s < 5.5 || phi_s > 0.6
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2008 recommended range',phi_s*100,char(37))        
    %        end              
    %    case 4
    %      if phi_s < 0.2 || phi_s > 0.8
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2008 allowable range',phi_s*100,char(37))        
    %      end                
    %end    

else     %% wet coil operation
    
    % water vapour saturation pressure at coil temperature, pascals
    Psat_coil = 611.21*exp((18.678 - Tcoil/234.5).*Tcoil/(257.14 + Tcoil));  
    
    % volume mixing ratio at coil temperature
    X = Psat_coil/P;

    % specific humidity at coil temperature, kg H20/kg humid air 
    Hcoil = X*M_H2O/(X*M_H2O + (1 - X)*M_dry);
    
    % specific humidity at supply temperature, kg H20/kg humid air
    Hs = H - (H - Hcoil)*(Tr - Ts)/(Tr - Tcoil);
    
    % specific enthalpy at supply temperature, joule/kg
    hs = Hs*c_evap + Hs*Ts*c_vap + Ts*c_a;
    
    % net enthalpy change, joule/kg
    delta_h = hr - hs;
    
    % change in specific humidity, kg H20/kg humid air
    delta_H = H - Hs;
    
    % volume mixing ratio of supply air
    Xs = Hs*M_dry./(Hs*(M_dry - M_H2O) + M_H2O);
    
    % water vapour saturation pressure of supply air, joule/kg
    Psat = 611.21*exp((18.678 - Ts/234.5).*Ts/(257.14 + Ts));
    
    % relative humidity of supply air
    phi_s = Xs*P/Psat;
   
    % water vapour partial pressure at return temperature, pascals
    Pvap = phi_s*Psat;

    % dry air partial pressure, pascals
    Pdry = P - Pvap;
    
    % specific heat capacity of air as a function of water vapour content, Joule/(kg.C)
    c_ah = c_a + 1132*Pvap./Pdry;
    
    % density of humid air, kg/m^3
    rho = Pdry./(Rdry*(Ts + 273)) + Pvap./(Rvap*(Ts + 273));
    
    %%%> required airflow for removal of sensible IT heat load, m^3/s
    f = Q./(rho.*c_ah*(Tr - Ts)*e_coil);
    
    %%%> total heat (sensible + latent), watts
    Qload = f*rho*delta_h;
    
    %%%> latent heat
    Qlat = f*rho*c_evap*(H - Hs);
    Qtest = Qload - Qlat;
    
    % calculates required fan power
    Pcrac = Pmax*(f/fmax).^3;
    if Pcrac > Pmax
        error('Required CRAH power exceeds maximum rated value')
    end
    Pcrac = Pcrac/e_fan;
    Pdis = Pcrac*(1 - e_fan);
    
    % adds heat load + fan power draw 
    Qload = Qload + Pdis;
    
    % warning generation
    %switch mode
    %    case 1
    %        if phi_s < 0.4 || phi_s > 0.55
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2004 recommended range',phi_s*100,char(37))        
    %        end
    %    case 2
    %        if phi_s < 0.2 || phi_s > 0.8
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2004 allowable range',phi_s*100,char(37))        
    %        end         
    %    case 3
    %        % supply air dewpoint calculation, Celsius
    %        gamma = a*Ts/(b + Ts) + log(phi_s);
    %        Tdew_s = b*gamma/(a - gamma);
    %        if Tdew_s < 5.5 || phi_s > 0.6
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2008 recommended range',phi_s*100,char(37))        
    %        end              
    %    case 4
    %      if phi_s < 0.2 || phi_s > 0.8
    %            warning('Supply relative humidity is %.2f%s, outside of ASHRAE''s 2008 allowable range',phi_s*100,char(37))        
    %      end                
    %end    

end

return


