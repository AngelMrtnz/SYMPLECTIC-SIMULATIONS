% Order 2 explicit Ralston method

function [T, Y] = method_ralston(f, t, ci)
    % This function performs the explicit Ralston method (a second-order Runge-Kutta method)
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
        k1 = f(t(i), y(:, i));                       % Compute the derivative (slope) at the current point
        k2 = f(t(i) + h*3/4, y(:, i) + h*k1*3/4);    % Compute the derivative at the point 3/4 of the way through the interval
        y(:, i+1) = y(:, i) + h*(k1 + 2*k2)/3;       % Update the solution using a weighted average of the slopes
    end

    T = t'; % Transpose time vector to column vector
    Y = y'; % Transpose solution matrix to have each row correspond to a time point
end
