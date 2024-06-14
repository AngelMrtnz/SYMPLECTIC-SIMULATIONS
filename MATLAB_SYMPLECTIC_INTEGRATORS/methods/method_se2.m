% Order 1 explicit Symplectic Euler 2 method
% For separable Hamiltonians only!

function [T, Y] = method_se2(f, t, ci)
    % This function performs the explicit Symplectic Euler 2 method (a first-order
    % symplectic integrator) to solve ODEs for separable Hamiltonians.
    % Inputs:
    % f - function handle representing the ODE (dy/dt = f(t, y))
    % t - vector of time points where the solution is computed
    % ci - column vector of initial conditions
    % Outputs:
    % T - column vector of time points (same as input t, transposed)
    % Y - matrix of solutions; each row corresponds to the solution at a time point in T

    m = length(t);        % Number of time points
    h = t(2) - t(1);      % Step size (assumes uniform spacing of time points)

    [fil, col] = size(ci);    % Get size of initial conditions vector
    if (col > fil)
        ci = ci';             % Ensure initial conditions are a column vector
    end

    n = length(ci) / 2;       % Half the length of the initial conditions vector

    y = zeros(2 * n, m);      % Initialize the solution matrix
    y(:, 1) = ci;             % Set initial conditions

    for i = 1:m-1
        dy = f(t(i+1), y(:, i));                     % Compute the derivative (slope) at the next time point
        y(1:n, i+1) = y(1:n, i) + h * dy(1:n);       % Update the first half of the state vector
        dy = f(t(i), y(:, i+1));                     % Compute the derivative at the current time point using updated state
        y(n+1:2*n, i+1) = y(n+1:2*n, i) + h * dy(n+1:2*n); % Update the second half of the state vector
    end

    T = t'; % Transpose time vector to column vector
    Y = y'; % Transpose solution matrix to have each row correspond to a time point
end


