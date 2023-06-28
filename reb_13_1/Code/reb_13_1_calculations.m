function reb_13_1_calculations
%reb_13_1_calculations Performs reb_13_1 calculations
    % calculate the responses
    [t,CA,CB,CY,CZ,T,r] = response;

    % plot CA vs. t
    figure;
    plot(t,CA,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('t (h)','FontSize', 14)
    ylabel('CA (mol/L)','FontSize', 14)
    %saveas(gcf,"activity_11_residuals_T.png")
    %legend({'Perfect Model/Data','Actual Model/Data'},'Location',...
    %    'northwest','FontSize',14)

    figure;
    plot(t,r,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('t (h)','FontSize', 14)
    ylabel('r (mol/L/h)','FontSize', 14)

end