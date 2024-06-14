% Order 2 explicit Stormer Verlet 21 method
% For separable Hamiltonians only!!

function [T, Y] = method_sv21(f, t, ci)
    % This function performs the explicit Stormer Verlet 21 method (a second-order
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
        % Stormer-Verlet Step 2
        dy = f(t(i), y(:, i));                         % Compute the derivative at the current time
        qm = y(1:n, i) + (h/2) * dy(1:n);              % Update the position halfway using the derivative
        dy = f(0.5 * (t(i+1) + t(i)), [qm; y(n+1:2*n, i)]); % Compute the derivative at the midpoint using updated position
        pm = y(n+1:2*n, i) + (h/2) * dy(n+1:2*n);      % Update the momentum halfway using the derivative

        % Stormer-Verlet Step 1
        dy = f(0.5 * (t(i+1) + t(i)), [qm; pm]);      % Compute the derivative at the midpoint using updated position and momentum
        y(n+1:2*n, i+1) = pm + (h/2) * dy(n+1:2*n);    % Update the momentum using the derivative
        dy = f(t(i+1), [qm; y(n+1:2*n, i+1)]);         % Compute the derivative at the next time using updated position and momentum
        y(1:n, i+1) = qm + (h/2) * dy(1:n);            % Update the position using the derivative
    end

    T = t'; % Transpose time vector to column vector
    Y = y'; % Transpose solution matrix to have each row correspond to a time point
end
