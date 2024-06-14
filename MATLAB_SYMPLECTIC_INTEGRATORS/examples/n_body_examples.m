function [G, m, ci] = n_body_examples(case_number)
%%% Returns initial conditions for different n-body cases
%%%
%%% Inputs:
%%%   case_number: String indicating the desired n-body case ('case 1', 'case 2', or 'case 3')
%%%
%%% Outputs:
%%%   G: Gravitational constant
%%%   m: Masses of the bodies
%%%   ci: Initial positions and velocities of the bodies

    if case_number == 'case 1' % 3-body chaotic case 1
        G = 1;
        m = [1, 5, 5]; % Masses of the bodies
        % Initial positions
        ob10 = [5, 1, 0];
        ob20 = [0, 0, 0];
        ob30 = [-5, 0, 1];
        % Initial velocities
        p_ob10 = [0, -1, 0];
        p_ob20 = [0, 0, 1];
        p_ob30 = [1, 1, 0];
        ci = [ob10, ob20, ob30, p_ob10, p_ob20, p_ob30]; % Concatenate initial positions and velocities
    end
    if case_number == 'case 2' % 3-body periodic (triple circle)
        G = 1;
        m = [1, 1, 1]; % Masses of the bodies
        % Initial conditions
        ci = [-0.0347, 1.1856, 0, 0.2693, -1.0020, 0, -0.2328, -0.5978, 0, 0.2495, -0.1076, 0, 0.2059, -0.9396, 0, -0.4553, 1.0471, 0];
    end
    if case_number == 'case 3' % Outer solar system
        G = 2.95912208286e-4; % Gravitational constant for the solar system
        % Masses of the bodies (Sun, Jupiter, Saturn, Uranus, Neptune, Pluto)
        m = [1.00000597682, 9.54786104043e-4, 2.85583733151e-4, 4.37273164546e-5, 5.17759138449e-5, 1 / (1.3 * 10^8)];
        % Initial positions and velocities of the bodies
        initial_positions = [0, 0, 0, -3.5023653, -3.8169847, -1.5507963, 9.0755314, -3.0458353, -1.6483708, 8.3101420, -16.2901086, -7.2521278, 11.4707666, -25.7294829, -10.8169456, -15.5387357, -25.2225594, -3.1902382];
        initial_velocities = [0, 0, 0, m(2) * 0.00565429, m(2) * -0.00412490, m(2) * -0.00190589, m(3) * 0.00168318, m(3) * 0.00483525, m(3) * 0.00192462, m(4) * 0.00354178, m(4) * 0.00137102, m(4) * 0.00055029, m(5) * 0.00288930, m(5) * 0.00114527, m(5) * 0.00039677, m(6) * 0.00276725, m(6) * -0.00170702, m(6) * -0.00136504];
        ci = [initial_positions, initial_velocities]; % Concatenate initial positions and velocities
    end
end
