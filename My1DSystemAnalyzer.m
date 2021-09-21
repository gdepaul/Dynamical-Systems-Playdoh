function [] = My1DSystemAnalyzer()

    h = 0.01;
    dxdt = @(x) -sin(x);
    epsilon = 10^(-2);
    x_bounds = [-2*pi 2*pi];
    t_bounds = [-2*pi 2*pi];

    T = x_bounds(1):h:x_bounds(2);
    X = t_bounds(1):h:t_bounds(2);

    Y = arrayfun(dxdt, X);
    plot(X,Y)
    hold on
    A = (Y < epsilon) & (Y > -1*epsilon);
    scatter(X(A), Y(A), 'filled');
    grid on
    hl = xlabel('$x$');
    set(hl, 'Interpreter', 'latex');
    hl = ylabel('$\dot{x}$');
    set(hl, 'Interpreter', 'latex');
    set(gca,'FontSize',20);

    figure()
    grid off

    [T, X] = meshgrid(T,X);

    Z1 = ones(size(T));
    Z2 = arrayfun(dxdt, X);

    streamslice(T,X, Z1, Z2);
    axis([t_bounds(1) t_bounds(2) x_bounds(1) x_bounds(2)]);

    hl = xlabel('$t$');
    set(hl, 'Interpreter', 'latex');
    hl = ylabel('Family of curves $x(t)$');
    set(hl, 'Interpreter', 'latex');
    set(gca,'FontSize',20);
    
end