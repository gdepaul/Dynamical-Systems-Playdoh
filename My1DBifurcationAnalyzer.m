function [] = My1DBifurcationAnalyzer()

    h = 0.01;
    U_ = -10:h:10;
    X_ = -5:h:5;
    dxdt = @(u, x) x^3 - u*x;
    d = 10;

    [U, X] = meshgrid(U_,X_);
    Z = arrayfun(dxdt, U, X);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Create New Figure
    %%%%%%%%%%%%%%%%%%%%%%%%%
    figure()

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
    [M,c] = contour(U,X,Z1,v,':black');
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
    [M,c] = contour(U,X,Z2,v,'black');
    c.LineWidth = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Axis Preferences
    %%%%%%%%%%%%%%%%%%%%%%%%%
    hl = xlabel('$\mu$');
    set(hl, 'Interpreter', 'latex');
    hl = ylabel('$x^*$');
    set(gca,'FontSize',20);
    set(hl, 'Interpreter', 'latex');
    hl = title('$x^*$ vs $\mu$');
    set(hl, 'Interpreter', 'latex');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    
end