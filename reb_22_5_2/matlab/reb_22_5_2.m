function reb_22_5_2()
% Calculations for Reaction Engineering Basics Example 22.5.2
    % given and known constants
    VFRin = 10.0; % L /min
    VFRout = 12.2; % L /min
    Cbefore = 2.0; % mmol /L

    % read the experimental data
    data_table = readtable('../data.csv','VariableNamingRule'...
        ,'preserve');
    t = table2array(data_table(:,1)); % min
    Cout = table2array(data_table(:,2)); % mmol /L

    % plot Cout vs. t
    figure;
    plot(t,Cout,'ok','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Elapsed Time (min)','FontSize', 14)
    ylabel('Tracer Concentration (mM)','FontSize', 14)
    saveas(gcf,"response.png")

    % calculate F
    Ndata = length(t);
    F = nan(Ndata,1);
    for i=1:Ndata
        F(i) = (VFRin*Cbefore - VFRout*Cout(i))/VFRin/Cbefore;
    end

    % plot F vs. lambda (t)
    figure;
    plot(t,F,'ok','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Age (min)','FontSize', 14)
    ylabel('F','FontSize', 14)
    saveas(gcf,"cum_age_dist.png")
end