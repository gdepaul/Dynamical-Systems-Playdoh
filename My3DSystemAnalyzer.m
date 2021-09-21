function [] = My3DSystemAnalyzer() 

    % Specify System
    x_dot = @(x, y, z) -2*x + 1*y + 0*z;
    y_dot = @(x, y, z) 0*x + -1*y + -2*z;
    z_dot = @(x, y, z) 0*x + 2*y + -1*z;

%     x_dot = @(x, y, z) -10*(y - x);
%     y_dot = @(x, y, z) 28*x - y - x*z;
%     z_dot = @(x, y, z) -8/3*z + x*y;
    
    % Constants
    epsilon = 10^(-3); % Degree of allowable similarity between fixed points identified numerically
    h = 0.1; % Resolution for vector field
    r = 100; % Bounding box radius around each fixed point
    num_curves = 50; % Number of curves to be drawn
    curve_res = 200; % Resolution of curves
    curve_length = 1; % Length of curve (relative to curve's derivative)

    % Identify Fixed Points
    syms x y z
    assume(x, 'real');
    assume(y, 'real');
    assume(z, 'real');
    f = x_dot(x,y,z);
    g = y_dot(x,y,z);
    h = z_dot(x,y,z);
    eqn1 = f == 0;
    eqn2 = g == 0;
    eqn3 = h == 0;
    eqns = [eqn1, eqn2, eqn3];
    S = solve(eqns,[x y z]);
    
    for i = 1:length(S.x) 
        J = [diff(f,x) diff(f,y), diff(f,z); diff(g,x) diff(g,y), diff(g,z); diff(h,x) diff(h,y), diff(h,z)];
        J = subs(subs(subs(J, x, S.x(i)), y, S.y(i)),z,S.z(i)); 
        [V,D] = eig(J);
        eigenvalues = diag(D);
        latex(V(:,1))
        latex(V(:,2))
        latex(V(:,3))
        for j=1:num_curves
            a = [double(S.x(i)) double(S.y(i)) double(S.z(i))] + [r*(rand()-0.5) r*(rand()-0.5) r*(rand()-0.5)];
            f = @(t,x)[x_dot(x(1), x(2), x(3))
                y_dot(x(1), x(2), x(3))
                z_dot(x(1), x(2), x(3))];
            t=linspace(0,curve_length,curve_res);
            [t,x]=ode45(f,t,a);
            cm = jet(numel(t));
            hold on
            for k = 1:size(x,1)-1
                hl = plot3([x(k,1),x(k+1,1)], [x(k,2),x(k+1,2)], [x(k,3),x(k+1,3)]);
                set(hl, 'LineStyle','-', 'Color',[0 0 0]);
            end
            hold off
            grid on
            view(50, 30)
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Axis Preferences
    %%%%%%%%%%%%%%%%%%%%%%%%%
    hl = xlabel('$x$');
    set(hl, 'Interpreter', 'latex');
    hl = ylabel('$y$');
    set(hl, 'Interpreter', 'latex');
    hl = zlabel('$z$');
    set(hl, 'Interpreter', 'latex');
    set(gca,'FontSize',20);
    
end