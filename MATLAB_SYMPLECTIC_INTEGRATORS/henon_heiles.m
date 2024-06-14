function henon_heiles
    clc
    clear all
    close all

    % Add necessary paths and configure plotting
    addpath('methods');
    addpath('plotting');
    plot_config;

    % Parameters
    lambda = 1;
    ci = [0.3, 0, 0.1, 0.5];
    ti = 0;
    tf = 100;
    h = 0.05;

    % Perform simulations using different numerical integration methods

    % Method 1: SV Faster
    tic;
    [T1, Y1] = method_sv_faster(@(t,y) f1(t,y,lambda), @(t,y) f2(t,y,lambda), [ti:h:tf], ci);
    H1 = Hamiltonian(T1, Y1, lambda);
    H1_rel = H1 - H1(1);
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']);

    % Method 2: RK3
    tic;
    [T2, Y2] = method_rk3(@(t,y) f(t,y,lambda), [ti:h:tf], ci);
    H2 = Hamiltonian(T2, Y2, lambda);
    H2_rel = H2 - H2(1);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']);

    % Method 3: RK4
    tic;
    [T3, Y3] = method_rk4(@(t,y) f(t,y,lambda), [ti:h:tf], ci);
    H3 = Hamiltonian(T3, Y3, lambda);
    H3_rel = H3 - H3(1);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']);

    % Method 4: ODE45
    options = odeset('RelTol', 1e-12);
    tic;
    [T4, Y4] = ode45(@(t,y) f(t,y,lambda), [ti:h:tf], ci, options);
    H4 = Hamiltonian(T4, Y4, lambda);
    H4_rel = H4 - H4(1);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']);

    % Method 5: ODE89
    tic;
    [T5, Y5] = ode89(@(t,y) f(t,y,lambda), [ti:h:tf], ci, options);
    H5 = Hamiltonian(T5, Y5, lambda);
    H5_rel = H5 - H5(1);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% ENERGY PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create a folder for the results if it doesn't exist
    relativeFolderPath = 'results/henon_heiles';
    if ~exist(relativeFolderPath, 'dir')
        mkdir(relativeFolderPath);
    end

    % Plot energy difference
    fig = figure(1);
    plot(T1, H1_rel, T2, H2_rel, T3, H3_rel, ':', T4, H4_rel, '-.', T5, H5_rel, '--');
    legend('S-V 1', 'RK-3', 'RK-4', 'ode45', 'ode89', 'FontSize', 12, 'Location', 'southwest');
    title('Difference in energy');
    xlabel('time');
    ylabel('H-H(0)');
    saveas(fig, fullfile(relativeFolderPath, 'energy_test.pdf'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% TOTAL TRAJECTORIES PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot total trajectory
    fig = figure(2);
    plot(Y1(:, 1), Y1(:, 2));
    title('Total trajectory');
    xlabel('x');
    ylabel('y');
    axis([-1, 1, -1, 1]);
    axis equal;
    saveas(fig, fullfile(relativeFolderPath, 'trajectory_test.pdf'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% H-H in Movement Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot H-H in movement (create a video)
    fig = figure(3);
    plot(Y1(1, 1), Y1(1, 2));
    axis([-1, 1, -1, 1]);
    videoFileName = sprintf('video_test.mp4');
    videoPath = fullfile(relativeFolderPath, videoFileName);
    video = VideoWriter(videoPath, 'MPEG-4');
    video.FrameRate = 10;
    open(video);
    step =1;
    grid on;
    for k = 1:step:length(Y1(:, 1))
        figure(3);
        km = max(1, k - 50);
        a = Y1(km:k, 1);
        b = Y2(km:k, 2);
        titol = strcat('t = ', num2str(T1(k)));
        plot(a, b);
        title(titol);
        xlabel('x');
        ylabel('y');
        axis([-1, 1, -1, 1]);
        axis square;
        grid on;
        pause(0.01);
        frame = getframe(fig);
        writeVideo(video, frame);
    end
    close(video);
end

% Differential equations of the system
function dy = f1(t, y, lambda)
    
    dy = zeros(2, 1);
    dy(1) = y(3);
    dy(2) = y(4);
end

function dy = f2(t, y, lambda)

    dy = zeros(2, 1);
    dy(1) = -y(1) - lambda * 2 * y(1) * y(2);
    dy(2) = -y(2) - lambda * (y(1)^2 - y(2)^2);
end

function dy = f(t, y, lambda)

    dy = zeros(4, 1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1) - lambda * 2 * y(1) * y(2);
    dy(4) = -y(2) - lambda * (y(1)^2 - y(2)^2);
end

% Hamiltonian function
function H = Hamiltonian(t, y, lambda)
    
    H = (1/2) * (y(:, 3).^2 + y(:, 4).^2) + (1/2) * (y(:, 1).^2 + y(:, 2).^2) + ...
        lambda * ((y(:, 1).^2) .* y(:, 2) - (y(:, 2).^3) / 3);
end



