function [Y_DA, fA] = reb_13_3_response(rxn_times)
% performs the calculations for the 2022 CE 329 Evaluation 3 

    % Known, constant quantities
    Rgas = 0.082057*1000; % cc atm/mol/K
    Renergy = 1.987; % cal/mol/K
    k01 = 3.34E9; % mol/cc/min/atm^2
    k02 = 4.99E9;
    E1 = 20500; % cal/mol
    E2 = 21800;
    dH1 = -6300; % cal/mol
    dH2 = -6900;
    CpA = 7.4; % cal/mol/K
    CpB = 8.6;
    CpD = 10.7;
    CpZ = 5.2;
    CpU = 10.3;
    V = 2000; % cc
    PA0 = 1.; % atm
    PB0 = 2.;
    T0 = 25 + 273.15; % K
    U = 0.6; % cal/cm^2/min/K
    A = 600; % cm^2
    Te = 30 + 273.15;

    % Calculated constant quantities
    nA0 = PA0*V/Rgas/T0;
    nB0 = PB0*V/Rgas/T0;
    P0 = PA0 + PB0;

    % Adjusted values
    nAdj = length(rxn_times);

    % Storage for responses
    Y_DA = nan(nAdj,1);
    fA = nan(nAdj,1);

    % Loop through the adjusted values
    for iTime = 1:nAdj
        % initial and final values
        t0 = 0.0;
        y0 = [
            nA0
            nB0
            0.0
            0.0
            0.0
            T0
            P0
        ];
        fVar = 0;
        fVal = rxn_times(iTime);

        % solver settings
        odesAreStiff = true;
        useMassMatrix = true;

        % solve the IVODEs
        [~, y, flag] = solveIVODEs(t0,y0,fVar,fVal,@desEqns,...
            odesAreStiff,useMassMatrix, @massMatrix);
        % check whether a solution was obtained
        if flag <= 0
            disp('The solution to the ODEs may not be accurate.')
        end

        % calculate the responses
        Y_DA(iTime) = y(end,3)/nA0;
        fA(iTime) = (nA0 - y(end,1))/nA0;
    end

    % function that evaluates the ODEs
    function dydx = desEqns(~, yVal)
        nA = yVal(1);
        nB = yVal(2);
        nD = yVal(3);
        T = yVal(6);
        PA = nA*Rgas*T/V;
        PB = nB*Rgas*T/V;
        PD = nD*Rgas*T/V;
        r1 = k01*exp(-E1/Renergy/T)*PA*PB;
        r2 = k02*exp(-E2/Renergy/T)*PD*PB;
        dydx = [
            -V*r1
            -V*(r1 + r2)
            V*(r1 - r2)
            V*(r1 + r2)
            V*r2
            U*A*(Te-T) - V*(r1*dH1 + r2*dH2)
            0.0
        ];
    end

    % function that evaluates the mass matrix
    function MM = massMatrix(~,yVal)
        Tcoeff = yVal(1)*CpA + yVal(2)*CpB + yVal(3)*CpD + yVal(4)*CpZ...
            + yVal(5)* CpU;
        minusRT = -Renergy*yVal(6);
        MM = [1, 0, 0, 0, 0, 0, 0
              0, 1, 0, 0, 0, 0, 0
              0, 0, 1, 0, 0, 0, 0
              0, 0, 0, 1, 0, 0, 0 
              0, 0, 0, 0, 1, 0, 0
              0, 0, 0, 0, 0, Tcoeff, -V*Renergy/Rgas
              minusRT, minusRT, minusRT, minusRT, minusRT, ...
              -Renergy*sum(yVal(1:5)), -V*Renergy/Rgas];
    end
end
