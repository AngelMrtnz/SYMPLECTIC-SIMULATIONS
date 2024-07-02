function oscillator
    % Clear workspace, command window, and close all figures
    clc
    clear all
    close all
    
    % Add path to the folder containing methods
    addpath('methods');
    addpath('examples');
    addpath('plotting');
    
    % Define system parameters
    m = 1;  % Mass
    k = 1;  % Spring constant
    q0 = 1; % Initial displacement
    p0 = 0; % Initial momentum
    ci = [q0,p0]; % Initial conditions
    ti = 0;      % Initial time
    tf = 100;    % Final time
    h = .01;     % Step size for numerical integration

    % Perform simulations using different numerical integration methods

    % Method 1: SV Faster
    tic;
    [T1, Y1] = method_sv_faster(@(t,y) f1(t,y,m,k), @(t,y) f2(t,y,m,k), [ti:h:tf], ci);
    H1 = Hamiltonian(T1, Y1, m,k);
    H1_rel = H1 - H1(1); % Relative change in Hamiltonian
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']);

    % Method 2: RK3
    tic;
    [T2, Y2] = method_rk3(@(t,y) f(t,y,m,k), [ti:h:tf], ci);
    H2 = Hamiltonian(T2, Y2, m,k);
    H2_rel = H2 - H2(1);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']);

    % Method 3: RK4
    tic;
    [T3, Y3] = method_rk4(@(t,y) f(t,y,m,k), [ti:h:tf], ci);
    H3 = Hamiltonian(T3, Y3,m,k);
    H3_rel = H3 - H3(1);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']);

    % Method 4: ODE45
    options = odeset('RelTol', 1e-12);
    tic;
    [T4, Y4] = ode45(@(t,y) f(t,y,m,k), [ti:h:tf], ci, options);
    H4 = Hamiltonian(T4, Y4,m,k);
    H4_rel = H4 - H4(1);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']);

    % Method 5: ODE89
    tic;
    [T5, Y5] = ode89(@(t,y) f(t,y,m,k), [ti:h:tf], ci, options);
    H5 = Hamiltonian(T5, Y5, m,k);
    H5_rel = H5 - H5(1);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% ENERGY PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create a folder for the results if it doesn't exist
    relativeFolderPath = 'results/oscillator';
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
    plot(Y1(:, 1),zeros(size(Y1)));
    title('Total trajectory');
    xlabel('x');
    ylabel('y');
    axis([-1, 1, -1, 1]);
    axis equal;
    saveas(fig, fullfile(relativeFolderPath, 'trajectory_test.pdf'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% oscillator in Movement Figure %%%%%%%%%%%%%%%%%%%%%%%

    % Plot oscillator in movement (create a video)
    fig = figure(3);
    plot(Y1(1, 1), 0);
    axis([-1, 1, -1, 1]);
    videoFileName = sprintf('video_test.mp4');
    videoPath = fullfile(relativeFolderPath, videoFileName);
    video = VideoWriter(videoPath, 'MPEG-4');
    video.FrameRate = 10;
    open(video);
    step = 5;
    grid on;
    for k = 1:step:length(Y1(:, 1))
        figure(3);
        km = max(1, k - 50);
        a = Y1(km:k, 1);
        titol = strcat('t = ', num2str(T1(k)));
        plot(a,zeros(size(a)));
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

% Define differential equations for numerical integration
function dy = f1(t,y,m,k)
    dy = zeros(1,1);
    dy(1) = y(2)/m; 
end

function dy = f2(t,y,m,k)
    dy = zeros(1,1);
    dy(1) = -k^2*m*y(1); 
end

function dy = f(t,y,m,k)
    dy = zeros(2,1);
    dy(1) = y(2)/m;        
    dy(2) = -k^2*m*y(1);   
end

% Calculate Hamiltonian of the system
function H = Hamiltonian(t, y, m, k)
    H = (0.5*y(:,2).^2)/m + 0.5*k^2*y(:,1).^2; % Kinetic + potential energy
end

