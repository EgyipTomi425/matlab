clc;
clear;
warning('off')

pkg load symbolic

function [H, df_dx, stationary_points] = hesse2d(f, x1, x2)
    % f példa: f_test = @(x1, x2) x1^3 + x2^3 - 3*x1 - 3*x2;

    % Függvény definiálása
    syms sx1 sx2
    f_sym = f(sx1, sx2);

    % Parciális deriváltak számítása
    df_dx1 = diff(f_sym, sx1);
    df_dx2 = diff(f_sym, sx2);

    % Másodrendű parciális deriváltak számítása
    d2f_dx1dx1 = diff(df_dx1, sx1);
    d2f_dx1dx2 = diff(df_dx1, sx2);
    d2f_dx2dx1 = diff(df_dx2, sx1);
    d2f_dx2dx2 = diff(df_dx2, sx2);

    % Hesse-mátrix összeállítása
    H = [d2f_dx1dx1, d2f_dx1dx2; d2f_dx2dx1, d2f_dx2dx2];

    % Alsó parciális derivált mátrix összeállítása
    df_dx = [df_dx1; df_dx2];

    % Stacionárius pontok kiszámítása
    stationary_points_x1 = solve(df_dx1 == 0, sx1);
    stationary_points_x2 = solve(df_dx2 == 0, sx2);

    % Descartes-szorzat létrehozása
    [Stationary_points_x1, Stationary_points_x2] = meshgrid(stationary_points_x1, stationary_points_x2);
    stationary_points = [Stationary_points_x1(:), Stationary_points_x2(:)];

    % Stacionárius pontokhoz tartozó maximum, minimum és nyeregpont értékek kiszámítása
    stationary_points_states = cell(size(stationary_points, 1), 1);
    for i = 1:size(stationary_points, 1)
        point = stationary_points(i, :);
        H_point = subs(H, [sx1, sx2], point);
        eigenvalues = eig(double(H_point));
        if all(eigenvalues < 0)
            stationary_points_states{i} = 'maximum';
        elseif all(eigenvalues > 0)
            stationary_points_states{i} = 'minimum';
        else
            stationary_points_states{i} = 'nyeregpont';
        end
    end

    % Stacionárius pontokhoz tartozó függvényértékek értékek kiszámítása
    stationary_points_values = zeros(size(stationary_points, 1), 1);
    for i = 1:size(stationary_points, 1)
        point = stationary_points(i, :);
        stationary_points_values(i) = double(subs(f_sym, [sx1, sx2], point));
    end

    % Stacionárius pontok típusának és értékének hozzáadása a stationary_points mátrixhoz
    stationary_points = [stationary_points, stationary_points_states, stationary_points_values];


    % Kiíratás
    disp("A két változós polinomiális függvény:\n");
    disp(f);

    disp("\nHesse-mátrix a tesztfüggvényre:\n");
    disp(H);

    disp("\nAlsó parciális derivált mátrix a tesztfüggvényre:\n");
    disp(df_dx);

    disp("\nStacionárius pontok a tesztfüggvényre:\n");
    disp(stationary_points);
    disp("\n\n");
end

% Tesztfüggvény definiálása
f_test = @(x1, x2) x1^3 + x2^3 - 3*x1 - 3*x2;

% Hesse-mátrix, alsó parciális derivált mátrix és stacionárius pontok számítása a tesztfüggvényre
[H_test, df_dx_test, stationary_points_test] = hesse2d(f_test, 'x1', 'x2');

% Feladatmegoldások a 8. feladatsorhoz:
f_test = @(x1, x2) 2*x1^3 - 3*x1^2 - 12*x1 + 1/3 * x2^3 - 4*x2;
[H_test, df_dx_test, stationary_points_test] = hesse2d(f_test, 'x1', 'x2');

f_test = @(x1, x2) 2*x1^3 + 3*x1^2-12*x1+2*x2^3-3*x2^2-12*x2;
[H_test, df_dx_test, stationary_points_test] = hesse2d(f_test, 'x1', 'x2');
