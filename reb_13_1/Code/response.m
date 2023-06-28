function [t,CA,CB,CY,CZ,T,r] = response
%response Response function for RE Basics Example 13-1
    % Given and known
    dH = -101.2E3; % J/mol
    k0 = 5.11e4 * 3600; % L/mol/h
    E = 74.8e3; % J/mol
    T0 = 180 + 273.15; % K
    V = 1900.0; % L
    CA0 = 2.9; % mol/L
    CB0 = 3.2; % mol/L
    Cp = 1.23 * 4.184; % J/g/K
    rho = 1.02 * 1000.0; % g/L
    tRxn = 2.0; % h
    R = 8.314; % J/mol/K

    ind0 = linspace(0,tRxn,100);
    dep0 = [CA0*V; CB0*V; 0.0; 0.0; T0];
    [t, depVar] = ode45(@desEqns, ind0, dep0);
    CA = depVar(:,1)/V;
    CB = depVar(:,2)/V;
    CY = depVar(:,3)/V;
    CZ = depVar(:,4)/V;
    T = depVar(:,5);
    r = k0*exp(-E/R./T).*CA.*CB;

    function ddt = desEqns(~,dep)
        nAlocal = dep(1);
        nBlocal = dep(2);
        Tlocal = dep(5);
        k = k0*exp(-E/R/Tlocal);
        rate = k*nAlocal/V*nBlocal/V;
        dAdt = -V*rate;
        dBdt = -V*rate;
        dYdt = V*rate;
        dZdt = V*rate;
        dTdt = -rate*dH/rho/Cp;
        ddt = [dAdt; dBdt; dYdt; dZdt; dTdt];
    end
end