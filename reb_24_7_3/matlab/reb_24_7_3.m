function reb_24_7_3()
% Calculations Reaction Engineering Basics Example 24.7.3
% revised 8/26/24
    % given and known constants
    Tfeed = 150 + 273.15; % K
    CAfeed = 2.0; % mol/L
    V = 500.0; % L
    VFR = 250.0; % L/h
    Tex = 180 + 273.15; % K
    A = 2.0; % m2
    U = 500E3; % cal/m2/h/K
    k0 = 1.14E9; % L/mol/h
    E = 16200.; % cal/mol
    Cp = 1.17E3; % cal/L/K
    dH = 18.2E3; % cal/mol
    R = 1.987; % cal/mol/K
    Dax = 3.0E-1;
    lambda = 1.0E1;

    % convert L to m^3
    liters_per_cu_meter = 1000.0;
    CAfeed = CAfeed*liters_per_cu_meter;
    V = V/liters_per_cu_meter;
    VFR = VFR/liters_per_cu_meter;
    k0 = k0/liters_per_cu_meter;
    Cp = Cp*liters_per_cu_meter;

    % calculated constants
    D = 4*V/A; % m
    L = A/pi/D; % m

    % CSTR residuals function
    function epsilon = CSTR_residuals(testSoln)
        nA1 = testSoln(1);
        nZ1 = testSoln(2);
        T1 = testSoln(3);
        
        % calculate unknowns
        nA0 = CAfeed*VFR;
        k = k0*exp(-E/R/T1);
        CA = nA1/VFR;
        r = k*CA^2;
        Q = U*A*(Tex - T1);

        % calculate the residuals
        epsilon_1 = nA0 - nA1 - r*V;
        epsilon_2 = -nZ1 + r*V;
        epsilon_3 = VFR*Cp*(T1 - Tfeed) + r*V*dH - Q;
        epsilon = [epsilon_1; epsilon_2; epsilon_3];
    end

    % CSTR model function
    function [nA1, nZ1, T1] = CSTR_unknowns()
        % guess
        guess = [
            VFR*CAfeed/2 % nA in mol/s
            VFR*CAfeed/2 % nB in mol/s
            Tex - 10.0 % T in K
        ];
        
        % solve the design equations
        [soln, flag, message] = solve_ates(@CSTR_residuals, guess);
        if flag <= 0
            disp(' ')
            disp('The design equations were NOT solved.')
            disp(['    Matlab provided the following message: '...
                ,message])
        end
        
        % extract results
        nA1 = soln(1);
        nZ1 = soln(2);
        T1 = soln(3);
    end

    % axial dispersion derivatives function
    function ddz = ad_derivatives(~,dep)
        % extract the dependent variables that are needed
        CA = dep(1);
        T = dep(3);
        wA = dep(4);
        wZ = dep(5);
        wT = dep(6);

        % calculate any other unknown quantities
        k = k0*exp(-E/R/T);
        r = k*CA^2;

        % evaluate the derivatives
        dCAdz = wA;
        dCZdz = wZ;
        dTdz = wT;
        dwAdz = -1/Dax*(-4*VFR*wA/pi/D^2 - r);
        dwZdz = -1/Dax*(-4*VFR*wZ/pi/D^2 + r);
        dwTdz = -1/lambda*(-4*VFR*Cp*wT/pi/D^2 + 4*U*(Tex - T)/D ...
            - r*dH);

        % combine the derivatives in a vector and return
        ddz = [dCAdz; dCZdz; dTdz; dwAdz; dwZdz; dwTdz];
    end

    % boundary values residuals function
    function epsilon = bv_residuals(depLow, depHi)
        % extract the dependent variables that are needed
        CA0 = depLow(1);
        CZ0 = depLow(2);
        T0 = depLow(3);
        wA0 = depLow(4);
        wZ0 = depLow(5);
        wT0 = depLow(6);
        wA1 = depHi(4);
        wZ1 = depHi(5);
        wT1 = depHi(6);

        % evaluate the residuals
        epsilon_1 = CAfeed + Dax*pi*D^2*wA0/4/VFR - CA0;
        epsilon_2 = Dax*pi*D^2*wZ0/4/VFR - CZ0;
        epsilon_3 = Tfeed + lambda*pi*D^2*wT0/4/VFR/Cp - T0;
        epsilon_4 = wA1;
        epsilon_5 = wZ1;
        epsilon_6 = wT1;

        % combine the residuals in a vector and return
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4 ...
            ; epsilon_5; epsilon_6];
    end

    % axial dispersion reactor function
    function [z, CA, CZ, T, wA, wZ, wT] = ad_profiles(guess)
        % solve the bvodes
        % set the upper and lower bounds of the independent variable
        xl = 0.0;
        xu = L;
        % solve the bvodes
        [z, dep] = solve_bvodes(xl, xu, guess, @ad_derivatives...
            , @bv_residuals);

        % extract the profiles
        CA = dep(1,:);
        CZ = dep(2,:);
        T = dep(3,:);
        wA = dep(4,:);
        wZ = dep(5,:);
        wT = dep(6,1);
    end

    % function that performs the calculations
    function perform_the_calculations()
        % solve the CSTR design equations
        [nA1, nZ1, T1] = CSTR_unknowns();

        % guess the axial dispersion profiles
        CAguess = nA1/VFR;
        CZguess = nZ1/VFR;
        Tguess = T1;
        wAguess = (nA1/VFR - CAfeed)/L;
        wZguess = -wAguess;
        wTguess = (T1 - Tfeed)/L;
        guess = [CAguess; CZguess; Tguess; wAguess; wZguess; wTguess];

        % solve the axial dispersion design equations
        [z, CA, ~, T, ~, ~, ~] = ad_profiles(guess);
        
        % display the results
        fA_cstr = 100.0*(CAfeed*VFR - nA1)/CAfeed/VFR;
        fA_ad = 100.0*(CAfeed - CA(end))/CAfeed;
        T_cstr = T1 - 273.15;
        T_ad = T(end) - 273.15;
        item = ["Axial Dispersion Diameter"...
            ; "Axial Dispersion Length"; "CSTR Conversion"...
            ; "Axial Dispersion Conversion"; "CSTR Temperature"...
            ; "Axial Dispersion Temperature"];
        value = [D; L; fA_cstr; fA_ad; T_cstr; T_ad];
        units = ["m"; "m"; "%"; "%"; "°C"; "°C"];
        resultsTable = table(item, value, units);
        disp(resultsTable)

        % plot the A concentration and temperature profiles
        fa = figure;
        hold on
        plot(z,CA/liters_per_cu_meter,'linewidth',2)
        yline(CAguess/liters_per_cu_meter, 'linewidth',2)
        hold off
        set(gca,'fontsize',14)
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Concentration of A (M)','FontSize', 14)
        saveas(gcf,"CA_profile.png")
        set(fa,'Units','Inches');
        pos = get(fa,'Position');
        set(fa,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fa,'CA_profile.pdf','-dpdf','-r0')

        fb = figure;
        hold on
        plot(z,T-273.15,'linewidth',2)
        yline(T1 - 273.15, 'linewidth',2)
        hold off
        set(gca,'fontsize',14)
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        saveas(gcf,"T_profile.png")
        set(fb,'Units','Inches');
        pos = get(fa,'Position');
        set(fb,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fb,'T_profile.pdf','-dpdf','-r0')

        % save the results
        writetable(resultsTable,'results.csv');
    end
    
    % perform the calculations
    perform_the_calculations();
end