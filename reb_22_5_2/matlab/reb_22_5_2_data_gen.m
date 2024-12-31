function reb_22_5_2_data_gen
% Generation of data for Reaction Engineering Basics Example 22.5.2
    % given
    VFRin = 10.0; % L /min
    VFRout = 12.2; % L /min
    Cbefore = 2.0; % mmol /L

    V = 5.0; % L

    nData = 41;
    t = nan(nData,1);
    Cout = nan(nData,1);
    tNext = 0.0;
    for i = 1:nData
        t(i) = tNext;
        tNext = tNext + 0.2;
        Cout(i) = VFRin*Cbefore/VFRout *exp(-VFRout*t(i)/V);
        % add noise
        Cout(i) = Cout(i) + (2*random('Normal',1,5) -1.0)*0.02/VFRout;
        Cout(i) = round(Cout(i),3,'significant');
        if Cout(i) < 0
            Cout(i) = 0.0;
        end
    end

    % plot Cout vs. t
    figure;
    plot(t,Cout,'ok','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Elapsed Time (min)','FontSize', 14)
    ylabel('Tracer Concentration (mM)','FontSize', 14)

    % save the data
    data_table = table(t,Cout);
    data_table.Properties.VariableNames = ["Elapsed Time (min)"...
        ,"Tracer Concentration (mM)"];
    disp(data_table)
    writetable(data_table,'../data.csv')
end