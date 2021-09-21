function [] = BasinsOfAttraction()

    n = 5;
    for i = 1:n
        e = 0.3*[cos(i*2*pi/n) sin(i*2*pi/n)];
        y0 = [1 0] + e;
        tspan = [0 200];
        sol = ode45(@vdp1,tspan,y0)
        x = linspace(0,200,25000);
        y = deval(sol,x);
        hold on
        plot(y(1,:), y(2,:))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Axis Preferences
        %%%%%%%%%%%%%%%%%%%%%%%%%
        hl = xlabel('$x$');
        set(hl, 'Interpreter', 'latex');
        hl = ylabel('$y$');
        set(gca,'FontSize',20);
        set(hl, 'Interpreter', 'latex');
        %hl = title('$x^*$ vs $h$');
        %set(hl, 'Interpreter', 'latex');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end

end

function dydt = vdp1(t,y)
    b = 0.01;
    dydt = -[y(2); -b*y(2) + y(1) - y(1)^3];
end