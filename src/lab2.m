function lab2(X, gamma)
    % Проверки на защите
    % X = [X X];
    % X = [X 100 100]; 
    n = length(X);

    % Выборочные оценки
    [mu, s2] = calc_select_params(X);
    
    fprintf("Выборочное среднее = %.3f\n", mu);
    fprintf("Исправленная выборочная дисперсия = %.3f\n", s2);

    % Вычисление доверительных интервалов

    alpha = (1 - gamma) / 2;

    [lower_m, upper_m] = calc_m_whithout_sigma_confint(X, alpha);

    fprintf("\ngamma-доверительный интервал для mu: (%.4f, %.4f)\n", lower_m, upper_m);

    [lower_sigma, upper_sigma] = calc_sigma_confint(X, alpha);

    fprintf("\ngamma-доверительный интервал для sigma: (%.4f, %.4f)\n", lower_sigma, upper_sigma);

    % Построение графиков

    mu_arr = zeros(n, 1);
    mu_line = zeros(n, 1);
    mu_lower = zeros(n, 1);
    mu_upper = zeros(n, 1);

    s2_arr = zeros(n, 1);
    s2_line = zeros(n, 1);
    sigma_lower = zeros(n, 1);
    sigma_upper = zeros(n, 1);

    mu_line(1 : n) = mu;
    s2_line(1 : n) = s2;
    for i = 1 : n
        X_ = X(1 : i);

        [mu_arr(i), s2_arr(i)] = calc_select_params(X_);
        [mu_lower(i), mu_upper(i)] = calc_m_whithout_sigma_confint(X_, alpha);
        [sigma_lower(i), sigma_upper(i)] = calc_sigma_confint(X_, alpha);
    end
    
    mu_plot(17, n, mu_line, mu_arr, mu_lower, mu_upper);
    figure();
    sigma_plot(17, n, s2_line, s2_arr, sigma_lower, sigma_upper);
end

function [mu, s2] = calc_select_params(X)
    n = length(X);
    
    mu = 0;
    s2 = 0;
    
    if (n > 0)
        mu = sum(X) / n;
    end
    if (n > 1)
        s2 = sum((X - mu) .^2) / (n - 1);
    end
end

% m - неизвестно, 
% sigma - неизвестно,
% Оценить m
function [l,u] = calc_m_whithout_sigma_confint(X, alpha)
    n = length(X);
    [mu, s2] = calc_select_params(X);
    q_st = tinv((1 - alpha), (n - 1));
    
    l = mu - (q_st * sqrt(s2) / sqrt(n));
    u = mu + (q_st * sqrt(s2) / sqrt(n));
end

% sigma - неизвестно
% Оценить sigma^2
function [l,u] = calc_sigma_confint(X, alpha)
    n = length(X);
    [~, s2] = calc_select_params(X);

    q_xi2_r = chi2inv((1 - alpha), (n - 1));
    q_xi2_l = chi2inv(alpha, (n - 1));
    
    l = s2 * (n - 1) / q_xi2_r;
    u = s2 * (n - 1) / q_xi2_l;
end

function mu_plot(startn, endn, mu_line, mu_arr, mu_lower, mu_upper)
    fprintf("Итоговый размах для mu = %.3f\n", mu_upper(end) - mu_lower(end))
    
    plot((startn : endn), mu_line(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), mu_arr(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), mu_upper(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), mu_lower(startn : endn), 'LineWidth', 1);
    hold on;
    
    grid on;
    xlabel("n");
    ylabel('\mu');

    legend('\mu\^(x_N)', '\mu\^(x_n)', '\mu^{-}(x_n)', '\mu_{-}(x_n)');
end

function sigma_plot(startn, endn, s2_line, s2_arr, sigma_lower, sigma_upper)
    fprintf("Итоговый размах для sigma = %.3f\n", sigma_upper(end) - sigma_lower(end));

    plot((startn : endn), s2_line(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), s2_arr(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), sigma_upper(startn : endn), 'LineWidth', 1);
    hold on;
    plot((startn : endn), sigma_lower(startn : endn), 'LineWidth', 1);
    hold on;
    
    grid on;
    xlabel("n");
    ylabel('\sigma');

    legend('S^2(x_N)', 'S^2(x_n)', '\sigma^{2 -}(x_n)', '\sigma^2_{-}(x_n)');
end
