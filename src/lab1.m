function lab1(X)
    X = sort(X);
    n = length(X);
    
    fprintf("а) Вычисление максимального значения Mmax " + ...
        "и минимального значения Mmin\n");
    
    Mmax = X(n);
    Mmin = X(1);
    
    fprintf("\nMmax = %.3f\n", Mmax);
    fprintf("Mmin  = %.3f\n", Mmin);
    
    fprintf("\nб) Вычисление размаха R\n");

    R = Mmax - Mmin;
    fprintf("\nR = %.3f\n", R);
    
    fprintf("\nв) Вычисление оценок Mu и S^2 " + ...
        "математического ожидания MX и дисперсии DX\n");
    
    Mu = sum(X) / n;
    S_square = sum((X - Mu) .^2) / (n - 1);
    
    fprintf("\nMu = %.3f\n", Mu);
    fprintf("S^2 = %.3f\n", S_square);
    
    fprintf("\nг) Группировка значений выборки в " + ...
        "m = [log2 n] + 2 интервала\n");
    
    m = floor(log2(n)) + 2; 
    fprintf("\nКол-во интервалов m = %3d:\n\n", m);
    
    delta = (X(n) - X(1)) / m;
    
    borders = Mmin : delta : Mmax;
    
    ni_arr = zeros(m, 1);
    
    cur_X_pos = 1;

    for i = 1 : m
        count = 0;
        lower_bound = borders(i);
        upper_bound = borders(i+1);
        is_last_interval = (i == m);
        
        while cur_X_pos <= length(X) && ( ...
                (is_last_interval && X(cur_X_pos) <= upper_bound) || ...
                (~is_last_interval && X(cur_X_pos) < upper_bound))
            if X(cur_X_pos) >= lower_bound
                count = count + 1;
            end
            cur_X_pos = cur_X_pos + 1;
        end
        
        if i == m
            fprintf(" %d. [%.3f; %.3f], кол-во элементов: %d\n", ...
                i, borders(i), borders(i+1), count);
        else
            fprintf(" %d. [%.3f; %.3f), кол-во элементов: %d\n", ...
                i, borders(i), borders(i+1), count);
        end
        
        ni_arr(i) = count;
    end

    fprintf("\nд) Построение гистограммы и " + ...
        "графика плотности нормального распределения\n\n"); 
    
    mid_intervals = zeros(m, 1);
    
    for i = 1 : m
        mid_intervals(i) = (borders(i) + borders(i + 1)) / 2;
    end
    
    column_values = zeros(m, 1);
    
    for i = 1 : m
        column_values(i) = ni_arr(i) / (n * delta);
    end
    
    figure;
    bar(mid_intervals, column_values, 1);
    hold on;
    
    x_coords = Mmin - 1 : 1e-3 : Mmax + 1;
    func_density_norm = normpdf(x_coords, Mu, sqrt(S_square)); 
    plot(x_coords, func_density_norm, 'LineWidth', 2);
    grid;
    legend('Гистограмма', 'Плотность N(\mu, \sigma^2)');
    
    fprintf("\nе) Построение эмпирической функции распределения " + ...
        "и нормальной ФР\n\n");   
    
    t_arr = [Mmin - 1, X, Mmax + 1];
    func_emperic = zeros(size(t_arr));
    
    for i = 1 : length(t_arr)
        func_emperic(i) = sum(X <= t_arr(i)) / n;
    end
    
    figure;
    stairs(t_arr, func_emperic, 'LineWidth', 1);
    hold on;
    
    func_norm = normcdf(x_coords, Mu, sqrt(S_square));
    plot(x_coords, func_norm, 'LineWidth', 1);
    grid;
    legend('Эмпирическая ФР', 'Нормальная ФР N(\mu, \sigma^2)');
end