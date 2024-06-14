% Order 2 explicit Heun method

function [T, Y] = method_heun(f, t, ci)
    % This function performs the explicit Heun method (also known as the improved Euler method)
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
        k1 = f(t(i), y(:, i));           % Compute the derivative (slope) at the current point
        k2 = f(t(i) + h, y(:, i) + h*k1); % Compute the derivative at the next point using the Euler step
        y(:, i+1) = y(:, i) + h*(k1 + k2)/2; % Update the solution using the Heun method
    end

    T = t'; % Transpose time vector to column vector
    Y = y'; % Transpose solution matrix to have each row correspond to a time point
end

