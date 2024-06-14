function pendulum
    % Clear workspace, command window, and close all figures
    clc
    clear all
    close all

    % Add path to the folder containing methods and examples
    addpath('methods');
    addpath('examples');

    % Load plot configuration settings
    plot_config;

    % Define system parameters
    m = 1;    % Mass
    l = 1;    % Length of the pendulum
    g = 9.8;  % Acceleration due to gravity
    q0 = pi/4; % Initial angle
    p0 = 0;    % Initial angular velocity
    ci = [q0,p0]; % Initial conditions
    ti = 0;    % Initial time
    tf = 5000; % Final time
    h = 0.01;  % Step size for numerical integration

    % Start the timer for function 1
    tic;
    [T1,Y1] = method_sv_faster(@(t,y) f1(t,y,m,l,g),@(t,y) f2(t,y,m,l,g), [ti:h:tf], ci);
    H1=Hamiltonian(T1,Y1,m,l,g);
    H1_rel=H1-H1(1); % Relative change in Hamiltonian
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']);

    % Start the timer for function 2
    tic;
    [T2,Y2] = method_rk3(@(t,y) f(t,y,m,l,g), [ti:h:tf], ci);
    H2=Hamiltonian(T2,Y2,m,l,g);
    H2_rel=H2-H2(1);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']);

    % Start the timer for function 3
    tic;
    [T3,Y3] = method_rk4(@(t,y) f(t,y,m,l,g), [ti:h:tf], ci);
    H3=Hamiltonian(T3,Y3,m,l,g);
    H3_rel=H3-H3(1);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']);  

    options = odeset('RelTol',1e-12);
    % Start the timer for function 4
    tic;
    [T4,Y4] = ode45(@(t,y) f(t,y,m,l,g), [ti:h:tf], ci,options);
    H4=Hamiltonian(T4,Y4,m,l,g);
    H4_rel=H4-H4(1);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']);

    % Start the timer for function 5
    tic;
    [T5,Y5] = ode89(@(t,y) f(t,y,m,l,g), [ti:h:tf], ci,options);
    H5=Hamiltonian(T5,Y5,m,l,g);
    H5_rel=H5-H5(1);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']);

    T = T1; % time for the video
    Y = Y1; % trajectory for the video
    x = x_c(l, Y(:,1));
	y = y_c(l, Y(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%Pendulum in movement figure%%%%%%%%%%%%%
    %{
    fig=figure(1);
    n = length(x);
    axis([-l-0.1,l+0.1,-0.1,2*l+0.1])

    % Specify the relative folder and filename
    relativeFolderPath = 'results/pendulum'; % Relative folder path (subfolder)
    videoFileName = sprintf('video_test.mp4');
    videoPath = fullfile(relativeFolderPath, videoFileName);
    video = VideoWriter(videoPath,'MPEG-4');

    % Set frame rate (optional)
    video.FrameRate = 10;
    % Open the video file
    open(video);
    step=10;
    grid on
    for k = 1:step:n
        figure(1)
        km = max(1, k-500);
        a = x(km:k);
        b = y(km:k);
        titol = strcat('t = ', num2str(T(k)));
        plot(a,b);
        hold on
        title(titol)
        line([0 x(k)], [l, y(k)], 'Color', 'black')
        hold off
        axis([-l-0.1,l+0.1,-0.1,2*l+0.1])
        axis square
        grid on
        pause(.01)
        % Capture the plot as a frame
        frame = getframe(fig);
        % Write the frame to the video file
        writeVideo(video, frame);
    end
    
    close(video);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% TOTAL TRAJECTORIES PLOT %%%%%%%%%%%%%%%%%%%%%%

    fig=figure(2); 
    plot(x,y)
    hold on
    title('Total trajectory')
    axis([-l-0.1,l+0.1,-0.1,2*l+0.1])
    xlabel('x');
    ylabel('y');
    axis equal;
    hold off

    % Specify the relative folder and filename
    relativeFolderPath = 'results/pendulum'; % Relative folder path (subfolder)
    fileName = 'trajectory_test.pdf';

    % Ensure the folder exists, if not, create it
    if ~exist(relativeFolderPath, 'dir')
        mkdir(relativeFolderPath);
    end

    % Full file path
    filePath = fullfile(relativeFolderPath, fileName);

    % Adjust the size of the paper to match the figure (optional)
    set(fig, 'Units', 'Inches');
    pos = get(fig, 'Position');
    set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

    % Save the figure as a PDF file
    saveas(fig, filePath);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%ENERGY PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    fig=figure(3);

    plot(T1,H1_rel,T2,H2_rel,T3,H3_rel,":", T4,H4_rel,"-.",T5,H5_rel,'--')
    legend('S-V 1','RK-3', 'RK-4','ode45','ode89', 'FontSize', 12, 'Location', 'southwest')
    title('Difference in energy')
    xlabel('time');
    ylabel('H-H(0)');
    hold off
    fileName = sprintf('energy_test.pdf');
    % Full file path
    filePath = fullfile(relativeFolderPath, fileName);

    % Save the figure as a SVG
    % Adjust the size of the paper to match the figure
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
    saveas(fig, filePath);

end

% Define differential equations for numerical integration
function dy = f(t,y,m,l,g)
    dy = zeros(2,1);
    dy(1) = y(2)/(m*l^2);             
    dy(2) = -m*g*l*sin(y(1));         
end

function dy1 = f1(t,y,m,l,g)
    dy1 =y(2)/(m*l^2);                
end

function dy2 = f2(t,y,m,l,g)
    dy2 = -m*g*l*sin(y(1));           
end

% Calculate Hamiltonian of the system
function H = Hamiltonian(t, y, m, l, g)
    H = (y(:,2).^2)/(2*m*l^2) + m*g*l*(1 - cos(y(:,1))); % Kinetic + potential energy
end

% Convert polar coordinates to Cartesian coordinates
function x=x_c(l,theta)
    x=l*sin(theta);
end

function y=y_c(l,theta)
    y=l*(1-cos(theta));
end
