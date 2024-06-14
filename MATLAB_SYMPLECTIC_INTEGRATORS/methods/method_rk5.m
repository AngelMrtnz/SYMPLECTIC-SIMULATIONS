% Order 5 explicit RK5 method

function [T, Y] = method_rk5(f, t, ci)
    % This function performs the explicit RK5 method (a fifth-order Runge-Kutta method)
    % to solve ODEs.
    % Inputs:
    % f - function handle representing the ODE (dy/dt = f(t, y))
    % t - vector of time points where the solution is computed
    % ci - column vector of initial conditions
    % Outputs:
    % T - column vector of time points (same as input t, transposed)
    % Y - matrix of solutions; each row corresponds to the solution at a time point in T

    n = length(t);        % Number of time points
    h = t(2) - t(1);      % Step size (assumes uniform spacing of time points)

    [fil, col] = size(ci);    % Get size of initial conditions vector
    if (col > fil)
        ci = ci';             % Ensure initial conditions are a column vector
    end

    y = zeros(length(ci), n); % Initialize the solution matrix
    y(:, 1) = ci;             % Set initial conditions

    for i = 1:n-1
        k1 = f(t(i), y(:, i));                           % Compute the derivative (slope) at the current point
        k2 = f(t(i) + h/4, y(:, i) + k1*h/4);            % Compute the derivative at 1/4 of the way through the interval
        k3 = f(t(i) + h/4, y(:, i) + h*(k1 + k2)/8);     % Compute the derivative at 1/4 of the way through using a combination of k1 and k2
        k4 = f(t(i) + h/2, y(:, i) + h*(-k2/2 + k3));    % Compute the derivative at the midpoint using a combination of k2 and k3
        k5 = f(t(i) + 3*h/4, y(:, i) + h*(3*k1/16 + 9*k4/16)); % Compute the derivative at 3/4 of the way through the interval
        k6 = f(t(i) + h, y(:, i) + h*(-3*k1/7 + 2*k2/7 + 12*k3/7 - 12*k4/7 + 8*k5/7)); % Compute the derivative at the end of the interval
        y(:, i+1) = y(:, i) + h*(7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)/90; % Update the solution using a weighted average of the slopes
    end

    T = t'; % Transpose time vector to column vector
    Y = y'; % Transpose solution matrix to have each row correspond to a time point
end
