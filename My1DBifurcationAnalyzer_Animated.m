function [] = My1DBifurcationAnalyzer_Animated()

    for a=-2:0.02:2
        h = 0.02;
        U_ = -2:h:2;
        X_ = -2:h:2;
        dxdt = @(u, x) x*(1-x) - u*x / (x + a);
        d = 2;

        [U, X] = meshgrid(U_,X_);
        Z = arrayfun(dxdt, U, X);

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Create New Figure
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %figure()
        clf

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw Unstable Region
        %%%%%%%%%%%%%%%%%%%%%%%%%
        v = [0,0];
        Z1 = Z;
        for i=1+d:length(X_)-d
            for j=1+d:length(U_)-d
                if(Z(i - d,j) > 0 && Z(i + d,j) < 0)
                    Z1(i,j) = NaN;
                end
                if(X(i) == -1) 
                    Z1(i,j) = NaN;
                end
            end
        end
        [M,c] = contourf(U,X,Z1,v,':black');
        c.LineWidth = 3;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw Stable Region
        %%%%%%%%%%%%%%%%%%%%%%%%%
        hold on
        v = [0,0];
        Z2 = Z;
        for i=1+d:length(X_)-d
            for j=1+d:length(U_)-d
                if(Z(i - d,j) < 0 && Z(i + d,j) > 0)
                    Z2(i,j) = NaN;
                end
            end
        end
        [M,c] = contourf(U,X,Z2,v,'black');
        c.LineWidth = 2;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Axis Preferences
        %%%%%%%%%%%%%%%%%%%%%%%%%
        hl = xlabel('$h$');
        set(hl, 'Interpreter', 'latex');
        hl = ylabel('$x^*$');
        set(gca,'FontSize',20);
        set(hl, 'Interpreter', 'latex');
        hl = title('$x^*$ vs $h$');
        set(hl, 'Interpreter', 'latex');
        %set(gca,'xtick',[])
        %set(gca,'ytick',[])

        pause(0.002)
    end

end