function reb_22_5_3()
% Calculations for Reaction Engineering Basics Example 22.5.3
    % conversion factors
    cm_per_m = 100.0;
    J_per_kJ = 1000.0;
    s_per_h = 3600.0;
    J_per_cal = 4.184;
    gPerCmS_per_cP = 0.01;
    atm_per_gPerSqSecPerCm = 9.869E-7;

    % given and known constants in cm, mol, J, s, K, atm, g
    yA0 = 0.75;
    yI0 = 0.25;
    P0 = 3.0; % atm
    T0 = 375 + 273.15; %K
    D = 2.5; % cm
    L = 8.0 * cm_per_m; % cm
    VFR0 = 100.0; % cm^3 /s
    mDot = 0.44; % g /s
    Tex = 375 + 273.15; % K
    U = 187.0 * J_per_kJ/s_per_h/cm_per_m^2; % J /s /cm^2 /K
    Dp = 0.25; % cm
    phi = 0.7;
    eps = 0.6;
    k0f = 9E17; % mol /cm^3 /s /atm
    Ef = 285.0 * J_per_kJ; % J /mol
    k0r = 4.09E-4; % mol /cm^3 /s /atm^4
    Er = 85.0 * J_per_kJ; % J /mol
    dH = 200.0 * J_per_kJ; % J /mol
    CpA = 11.7 * J_per_cal; % J /mol /K
    CpY = 8.3 * J_per_cal; % J /mol /K
    CpZ = 4.3 * J_per_cal; % J /mol /K
    CpI = 5.8 * J_per_cal; % J /mol /K
    mu = 0.27 * gPerCmS_per_cP; % g /cm /s
    Rpv = 82.057; % cm^3 atm /mol /K
    Ren = 8.314; % J /mol /K

    % calculated constants
    nA0 = yA0*P0*VFR0/Rpv/T0;
    nI0 = yI0*P0*VFR0/Rpv/T0;
    A = pi*D^2/4;
    G = mDot/A;

    % PFR derivatives function
    function ddz = PFR_derivatives(~,dep)
        % extract the dependent variables
        nA = dep(1);
        nY = dep(2);
        nZ = dep(3);
        nI = dep(4);
        T = dep(5);
        P = dep(6);

        % calculate the rate
        kf = k0f*exp(-Ef/Ren/T);
        kr = k0r*exp(-Er/Ren/T);
        PA = nA*P/(nA + nY + nZ + nI);
        PY = nY*P/(nA + nY + nZ + nI);
        PZ = nZ*P/(nA + nY + nZ + nI);
        r = kf*PA - kr*PY*PZ^3;

        % calculate the density
        rho = mDot*P/(nA + nY + nZ + nI)/Rpv/T;

        % evaluate the derivatives
        dnAdz = -A*r;
        dnYdz = A*r;
        dnZdz = 3*A*r;
        dnIdz = 0;
        dTdz = (pi*D*U*(Tex - T) - A*r*dH)/(nA*CpA + nY*CpY + nZ*CpZ...
            + nI*CpI);
        dPdz = -(1-eps)/eps^3*G^2/rho/phi/Dp*(150*(1-eps)*mu/phi/Dp/G ...
            + 1.75)*atm_per_gPerSqSecPerCm;

        % combine the derivatives in an array and return
        ddz = [dnAdz; dnYdz; dnZdz; dnIdz; dTdz; dPdz];
    end

    % PFR reactor function
    function [z, nA, nY, nZ, nI, T, P] = PFR_profiles()
        % initial values
        ind_0 = 0.0;
        dep_0 =[nA0; 0.0; 0.0; nI0; T0; P0];

        % stopping criterion
        f_var = 0;
        f_val = L;

        % solve the PFR design equations
        odes_are_stiff = false;
        [z, dep, flag, message] = solve_ivodes(ind_0, dep_0, f_var...
            , f_val, @PFR_derivatives, odes_are_stiff);
        if flag <= 0
            disp()
            disp(['Error solving ODEs: ', message])
        end

        % extract the profiles
        nA = dep(:,1);
        nY = dep(:,2);
        nZ = dep(:,3);
        nI = dep(:,4);
        T = dep(:,5);
        P = dep(:,6);
    end

    % quantities of interest
    function quantities_of_interest()
        % solve the PFR design equations
        [z, nA, nY, nZ, nI, T, P] = PFR_profiles();

        % calculate 1/VFR
        VFR = Rpv*(nA + nY + nZ + nI).*T./P;
        reciprocalVFR = 1./VFR;

        % calculate outlet conditions
        conversion = 100*(nA(1) - nA(end))/nA(1);
        Tout = T(end) - 273.15;

        % calculate t-bar
        t_bar = A*trapz(z,reciprocalVFR);
        t_bar_const_VFR = A*L/VFR0;

        % report and save the results
        item = ["t-bar with variable VFR";"t-bar with constant VFR"...
            ;"Outlet Conversion"; "Outlet Temperature"...
            ; "Outlet Pressure"];
        value = [t_bar; t_bar_const_VFR; conversion; Tout; P(end)];
        units = ["s"; "s"; "%"; "°C"; "atm"];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'results.csv')

        % plot VFR vs. z
        figure;
        plot(z,VFR,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (cm)','FontSize', 14)
        ylabel('Volumetric Flow Rate (cm^3 s^-^1)','FontSize', 14)
        saveas(gcf,"VFR_profile.png")
    end

    % perform the calculations
    quantities_of_interest();
end