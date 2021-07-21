%% Graph epidemic simulation produced by IH model
% Kate Royce
% October 2020

function produce_graph(species_name, t, x)
    figure;
    hold on;
          %S  (blue)          I  (red)           R (green)
    s1 = subplot(1,3,1);
    plot(t,x(:,1),'b*-',t,x(:,2),'r*-',t,x(:,3),'g*-'); %wild
    xlabel('t (days)'); ylabel('proportion'); title('Reservoir');
    set(s1, 'XLim', [0 120]);
    set(s1, 'YLim', [0 1]);
    set(gca,'FontSize',20);

    s2 = subplot(1,3,2);
    plot(t,x(:,4),'b*-',t,x(:,5),'r*-',t,x(:,6),'m*-',t,x(:,7),'g*-'); %domestic
    xlabel('t (days)'); ylabel('proportion'); title('Intermediate');
    legend({'susceptible','wildtype infected', 'transmissible infected', 'recovered'}, 'Location', 'north');
    set(gca,'FontSize',20);
    set(s2, 'XLim', [0 120]);
    set(s2, 'YLim', [0 1]);

    s3 = subplot(1,3,3);
    plot(t,x(:,8),'b*-',t,x(:,9),'m*-',t,x(:,10),'g*-'); %human
    xlabel('t (days)'); ylabel ('proportion'); title('Human');
    sgtitle(species_name,'FontSize',25);
    set(gca,'FontSize',20);
    set(s3, 'XLim', [0 120]);
    set(s3, 'YLim', [0 1]);
end