function [rNet, t, pctConv, CelsiusT] = reb_13_4_response(mEx)
% performs the calculations for the 2022 CE 329 Evaluation 3 

    % Known, constant quantities
    Renergy = 1.987; % cal/mol/K
    k01 = 2.59E9; % /min
    E1 = 16500; % cal/mol
    dH1 = -22200; % cal/mol
    CA0 = 2.0; % M
    T0 = 23 + 273.15; % K
    Cp = 440.0; % cal/L/K
    V = 4.0; % L
    Vex = 0.5; % L
    Aex = 0.6; % ft^2
    Uex = 1.13E4/60.0; % cal/ft^2/min/K
    TexIn = 20 + 273.15; % K
    rhoEx = 1.0; % g/cm^3
    CpEx = 1.0; % cal/g/K
    Ucoil = 3.8E4/60.0; % cal/ft^2/min/K
    Acoil = 0.23; % ft^2
    Tcoil = 120 + 273.15; % K
    Tex0 = 23 + 273.15; % K
    mEx0 = 0.0;
    T1 = 50 + 273.15; % K
    Tf = 25 + 273.15; % K
    tTurn = 25.0; % min

    % Calculated constant quantities
    nA0 = CA0*V;

    % set the phase
    heating = true;
    
    % initial and final values for heating phase
    t0 = 0.0;
    y0 = [
        nA0
        0.0
        T0
        Tex0
    ];
    fVar = 3;
    fVal = T1;

    % solver settings
    odesAreStiff = false;
    useMassMatrix = false;

    % solve the IVODEs
    [tHeat, yHeat, flag] = solveIVODEs(t0,y0,fVar,fVal,@desEqns,...
        odesAreStiff,useMassMatrix, 0);
    % check whether a solution was obtained
    if flag <= 0
        disp('The solution to the ODEs may not be accurate.')
    end
    
    % set the phase
    heating = false;

    % initial and final values for cooling phase
    t0 = tHeat(end);
    y0 = yHeat(end,:);
    fVal = Tf;

    % solve the IVODEs
    [t, y, flag] = solveIVODEs(t0,y0,fVar,fVal,@desEqns,...
        odesAreStiff,useMassMatrix, 0);
    % check whether a solution was obtained
    if flag <= 0
        disp('The solution to the ODEs may not be accurate.')
    end

    % concatenate the results
    t = [tHeat;t];
    y = [yHeat;y];

    % calculate the responses
    rNet = y(end,2)/(t(end) + tTurn);
    pctConv = 100.0*(nA0 - y(:,1))/nA0;
    CelsiusT = y(:,3) - 273.15;

    % function that evaluates the ODEs
    function ddt = desEqns(~, yVal)
        nA = yVal(1);
        T = yVal(3);
        Tex = yVal(4);
        k1 = k01*exp(-E1/Renergy/T);
        r1 = k1*nA/V;
        if heating
            mDotEx = mEx0;
            Qcoil = Ucoil*Acoil*(Tcoil - T);
        else
            mDotEx = mEx;
            Qcoil = 0;
        end
        Qex = Uex*Aex*(Tex - T);
        ddt = [
            -V*r1
            V*r1
            (Qex + Qcoil - V*r1*dH1)/(V*Cp)
            -((Qex + mDotEx*CpEx*(Tex - TexIn))/(rhoEx*Vex*CpEx))
        ];
    end
end
