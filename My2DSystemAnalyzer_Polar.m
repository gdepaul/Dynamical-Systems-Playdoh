function [] = My2DSystemAnalyzer_Polar() 
    
    % Symbolic System
	symbolic = 0;

    % Polar Coordinate System
    drdt = @(r,theta) r*(1 - r^2)*(4 - r^2);
    dthetadt = @(r,theta) 2 - r^2;
    
    % Specify System
    x_dot = @(x, y) drdt(sqrt(x^2 + y^2), angle([x,y])) * x / sqrt(x^2 + y^2) - y*dthetadt(sqrt(x^2 + y^2),angle([x,y])); 
    y_dot = @(x, y) drdt(sqrt(x^2 + y^2), angle([x,y])) * y / sqrt(x^2 + y^2) + x*dthetadt(sqrt(x^2 + y^2),angle([x,y]));
    
    % Constants
    epsilon = 10^(-3); % Degree of allowable similarity between fixed points identified numerically
    h = 0.01; % Resolution for vector field
    r = 2; % Bounding box radius around each fixed point
    vector_normalizer = 3; % Divide the vectors drawn for decoupled fixed points
    
    % Draw Real Vectos
    draw_real_eigenvectors = 1;
    
    % Draw Nullclines
    draw_nullclines = 0;
    
    % Override Plot Bounds
    override_bounds = 0;
    new_plot_bounds = [-2 2 -3 3];
    
    % Override Numerical Solver Params
    % NOTICE: Turn off when dealing with purely polynomial equations
    numerical_override = 0;
    numerical_search_range = [ -2 2; -2 2];
    number_of_possible_fixed_points = 11;

    % Identify Fixed Points
    syms x y x_0 y_0 t
    assume(x, 'real');
    assume(y, 'real');
    f = x_dot(x,y);
    g = y_dot(x,y);
    eqn1 = f == 0;
    eqn2 = g == 0;
    eqns = [eqn1, eqn2];
    if numerical_override == 1
        S = vpasolve(eqns, [x y], numerical_search_range);
        clear sol
        j = 1;
        for i = 1:(number_of_possible_fixed_points^2)
            temp = vpasolve(eqns, [x y], numerical_search_range, 'random', true);
            if ~isempty(temp.x)
                sol(i, 1) = temp.x;
                sol(i, 2) = temp.y;
                j = j + 1;
            end
        end
        
        % Sort Solution Array
        if j > 1
            sol = unique(sol, 'rows');

            % Identify Repeated Rows
            clear A
            A = zeros(length(sol));
            for i = length(sol):-1:2
                if abs(sol(i, 1) - sol(i-1, 1)) < epsilon && abs(sol(i, 2) - sol(i-1, 2)) < epsilon
                    A(i) = 1;
                end
            end

            % Isolate Unique Rows
            unique_x_values = zeros(length(sol) - length(A));
            unique_y_values = zeros(length(sol) - length(A));
            j = 1;
            for i = 1:length(sol)
                if A(i) == 0
                    unique_x_values(j) = sol(i,1);
                    unique_y_values(j) = sol(i,2);
                    j = j + 1;
                end
            end

            % Place In Struct
            S.x = unique_x_values;
            S.y = unique_y_values;
        end
    else
        S = solve(eqns,[x y]);
    end
    
    % Perform Calculations About Each Fixed Point
    if isempty(S.x)
        fprintf('No steady states found\n');
        min_X = -1;
        max_X = 1;
        min_Y = -1;
        max_Y = 1;
    else
        min_X = 100;
        max_X = -100;
        min_Y = 100;
        max_Y = -100;
        for i = 1:length(S.x)   
%             fprintf('Applying Linear Stability Analysis: \n');
%             
            if symbolic == 0
                if S.x(i) > max_X
                    max_X = S.x(i);
                end
                if S.x(i) < min_X
                    min_X = S.x(i);
                end
                if S.y(i) > max_Y
                    max_Y = S.y(i);
                end
                if S.y(i) < min_Y
                    min_Y = S.y(i);
                end
            end
            J = [diff(f,x) diff(f,y); diff(g,x) diff(g,y)];
            J = subs(subs(J, x, S.x(i)), y, S.y(i))
            [V,D] = eig(J);
            fprintf('Steady State:');
            disp([S.x(i) S.y(i)]);
            eigenvalues = diag(D)
            
            if symbolic ~= 1
                if eigenvalues == real(eigenvalues)

                    if eigenvalues(1) ~= eigenvalues(2)
                        % Distinct Real Eigenvalues
                        if min(eigenvalues) > 0
                            fprintf('(%.2f, %.2f) is an unstable node \n',S.x(i), S.y(i));
                        elseif max(eigenvalues) < 0
                            fprintf('(%.2f, %.2f) is a stable node \n',S.x(i), S.y(i));
                        else
                            fprintf('(%.2f, %.2f) is a saddle point \n',S.x(i), S.y(i));
                        end
                    else
                        eigenvectors_size = size(V);
                        % Repeated Eigenvalues
                        if eigenvectors_size(2) == 1
                            % Emits a Single Eigenvectors
                            if eigenvalues(1) > 0
                                fprintf('(%.2f, %.2f) is an unstable degenerate node \n',S.x(i), S.y(i));
                            else
                                fprintf('(%.2f, %.2f) is a stable degenerate node \n',S.x(i), S.y(i));
                            end
                        else
                            % Two Linearly Independent Eigenvectors
                            fprintf('(%.2f, %.2f) is a star node \n',S.x(i), S.y(i));
                        end
                    end

                    if draw_real_eigenvectors == 1
                        eigenvectors_size = size(V);

                        myColor = [rand() rand() rand()];
                        for k=1:eigenvectors_size(2)
                            [V(1, k) V(2,k)]
                            quiver(double(S.x(i)), double(S.y(i)), V(1, k)/vector_normalizer, V(2,k)/vector_normalizer,'linewidth',5,'color',myColor)
                            hold on
                            quiver(double(S.x(i)), double(S.y(i)), -V(1, k)/vector_normalizer, -V(2,k)/vector_normalizer,'linewidth',5,'color',myColor)
                            hold on
                        end
                    end
                else
                    % Complex Eigenvalues
                    real_component = real(eigenvalues(1));
                    if real_component > 0
                        fprintf('(%.2f, %.2f) is an unstable spiral \n',S.x(i), S.y(i));
                    elseif real_component < 0
                        fprintf('(%.2f, %.2f) is a stable spiral \n',S.x(i), S.y(i));
                    else
                        fprintf('(%.2f, %.2f) is a center \n',S.x(i), S.y(i));
                    end
                end
            end
        end
    end
    
    if symbolic ~= 1
    
        % Account for resolution
        min_X = min_X - r; 
        max_X = max_X + r;
        min_Y = min_Y - r; 
        max_Y = max_Y + r;

        % In case bounds are still bad
        if override_bounds == 1
            min_X = new_plot_bounds(1);
            max_X = new_plot_bounds(2);
            min_Y = new_plot_bounds(3);
            max_Y = new_plot_bounds(4);
        end

        % Draw Vector Field
        X_ = double(min_X):h:double(max_X);
        Y_ = double(min_Y):h:double(max_Y);
        [X, Y] = meshgrid(X_,Y_);
        Z_1 = arrayfun(x_dot, X, Y);
        Z_2 = arrayfun(y_dot, X, Y);
        streamslice(X,Y, Z_1, Z_2);
        axis([double(min_X) double(max_X) double(min_Y) double(max_Y)])

        if draw_nullclines == 1
            v = [0,0];
            hold on
            [M,c] = contour(X,Y,Z_1,v,':red');
            c.LineWidth = 3;
            hold on
            [M,c] = contour(X,Y,Z_2,v,':green');
            c.LineWidth = 3;
        end

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
        %set(gca,'xtick',[])
        %set(gca,'ytick',[])
    end
    
end