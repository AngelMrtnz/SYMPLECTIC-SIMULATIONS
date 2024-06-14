function two_body
    % Clear workspace, command window, and close all figures
    clc
    clear all
    close all
    
    % Add path to the folder containing methods and examples
    addpath('methods');
    addpath('plotting');
    addpath('examples');
    
    % Load plot configuration settings
    plot_config;

    % Define system parameters
    m1 = 1;   % Mass of first body
    m2 = 5;   % Mass of second body
    G = 1;    % Gravitational constant
    ci = [5,0,0,0,0,-1,0,1]; % Initial conditions [x1, y1, x2, y2, vx1, vy1, vx2, vy2]
    ti = 0;   % Initial time
    tf = 1000; % Final time
    h = 0.01; % Step size for numerical integration

    % Start the timer for function 1
    tic;
    [T1,Y1] = method_sv_faster(@(t,y) f1(t,y,m1,m2,G), @(t,y) f2(t,y,m1,m2,G), [ti:h:tf], ci);
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']); 
    H1=Hamiltonian(T1,Y1,m1,m2,G);
    L1=angular_momentum(T1,Y1,m1,m2,G);
    H1_rel=H1-H1(1);

    % Start the timer for function 2
    tic;
    [T2,Y2] = method_rk3(@(t,y) f(t,y,m1,m2,G), [ti:h:tf], ci);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']); 
    H2=Hamiltonian(T2,Y2,m1,m2,G);
    L2=angular_momentum(T2,Y2,m1,m2,G);
    H2_rel=H2-H2(1);

    % Start the timer for function 3
    tic;
    [T3,Y3] = method_rk4(@(t,y) f(t,y,m1,m2,G), [ti:h:tf], ci);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']); 
    H3=Hamiltonian(T3,Y3,m1,m2,G);
    L3=angular_momentum(T3,Y3,m1,m2,G);
    H3_rel=H3-H3(1);

    options = odeset('RelTol',1e-12);
    % Start the timer for function 4
    tic;
    [T4,Y4] = ode45(@(t,y) f(t,y,m1,m2,G), [ti:h:tf], ci,options);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']); 
    H4=Hamiltonian(T4,Y4,m1,m2,G);
    L4=angular_momentum(T4,Y4,m1,m2,G);
    H4_rel=H4-H4(1);

    % Start the timer for function 5
    tic;
    [T5,Y5] = ode89(@(t,y) f(t,y,m1,m2,G), [ti:h:tf], ci,options);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']); 
    H5=Hamiltonian(T5,Y5,m1,m2,G);
    L5=angular_momentum(T5,Y5,m1,m2,G);
    H5_rel=H5-H5(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% ENERGY AND ANGULAR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig=figure(3);
    set(fig, 'Position', [100, 100, 1000, 500]); % [left, bottom, width, height]

    subplot(1,2,1)
    plot(T1,H1_rel,T2,H2_rel,'--',T3,H3_rel,":", T4,H4_rel,"-.",T5,H5_rel)
    title('Difference in energy')
    xlabel('time');
    ylabel('H-H(0)');
    legend('S-V 1','RK-3', 'RK-4','ode45','ode89','FontSize', 10,'Location', 'southwest')

    subplot(1,2,2)
    plot(T1, abs(L1-L1(1)), T2, abs(L2-L2(1)),'--',T3, abs(L3-L3(1)),':', T4, abs(L4-L4(1)),'-.',T5, abs(L5-L5(1)))
    title('Difference angular momentum')
    xlabel('time');
    ylabel({'L-L(0) (norm)'});
    legend('S-V 1','RK-3', 'RK-4','ode45','ode89','FontSize', 10,'Location', 'northwest')

    % Adjust the position of the first subplot
    h1 = subplot(1, 2, 1);
    pos1 = get(h1, 'Position'); % Get the current position
    pos1(1)=pos1(1)-0.06;
    pos1(3) = pos1(3) * 1.1; % Increase the width
    set(h1, 'Position', pos1); % Set the new position
    
    % Adjust the position of the second subplot
    h2 = subplot(1, 2, 2);
    pos2 = get(h2, 'Position'); % Get the current position
    pos2(1) = pos2(1) + 0; % Shift the subplot to the left to avoid overlap
    pos2(3) = pos2(3) * 1.1; % Increase the width
    set(h2, 'Position', pos2); % Set the new position

    % Specify the relative folder and filename
    relativeFolderPath = 'results/two_body'; % Relative folder path (subfolder)
    fileName = 'energy_2_test.pdf';
    
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
    
    % Save the figure as an pdf file
    saveas(fig, filePath);
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax=plot_n_bodies_2d(T3,Y3,'two_body','2_rk5');
    plot_n_bodies_2d(T4,Y4,'two_body','2_midpoint',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T5,Y5,'two_body','2_ralston',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T2,Y2,'two_body','2_rk3',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T1,Y1,'two_body','2_sv21',1,10,.001,ax);

end

%Differential equations of the problem
function dy=f1(t,y,m1,m2,G)
    dy=zeros(4,1);
    dy(1)=y(5)/m1;
    dy(2)=y(6)/m1;
    dy(3)=y(7)/m2;
    dy(4)=y(8)/m2;
end

function dy=f2(t,y,m1,m2,G)
    dy=zeros(4,1);
    
    dif_x=y(3)-y(1);
    dif_y=y(4)-y(2);
    abs=(dif_x^2+dif_y^2)^(1/2);
    
    dy(1)= G*m1*m2*dif_x/(abs^3);
    dy(2)= G*m1*m2*dif_y/(abs^3);
    dy(3)=-G*m1*m2*dif_x/(abs^3);
    dy(4)=-G*m1*m2*dif_y/(abs^3);
end

function dy=f(t,y,m1,m2,G)
    dy=zeros(8,1);
    
    dif_x=y(3)-y(1);
    dif_y=y(4)-y(2);
    abs=(dif_x^2+dif_y^2)^(1/2);

    dy(1)=y(5)/m1;
    dy(2)=y(6)/m1;
    dy(3)=y(7)/m2;
    dy(4)=y(8)/m2;
    
    dy(5)= G*m1*m2*dif_x/(abs^3);
    dy(6)= G*m1*m2*dif_y/(abs^3);
    dy(7)=-G*m1*m2*dif_x/(abs^3);
    dy(8)=-G*m1*m2*dif_y/(abs^3);

end

%Calculates the Hamiltonian function
function H=Hamiltonian(t,y,m1,m2,G)
    dif_x=y(:,3)-y(:,1);
    dif_y=y(:,4)-y(:,2);
    abs=(dif_x.^2+dif_y.^2).^(1/2);

    H=(1/(2*m1))*(y(:,5).^2+y(:,6).^2)+(1/(2*m2))*(y(:,7).^2+y(:,8).^2)-G*m1*m2*(abs.^(-1));
end

function L = angular_momentum(t,y, m1, m2,G)
    % Calculate the angular momentum for the two-body system
    r1 = y(:,1:2);
    r2 = y(:,3:4);
    p1 = y(:,5:6);
    p2 = y(:,7:8);
    
    % Angular momentum in 2D
    L1 = r1(:,1).*p1(:,2) - r1(:,2).*p1(:,1); % L1z = r1 cross p1 in z direction
    L2 = r2(:,1).*p2(:,2) - r2(:,2).*p2(:,1); % L2z = r2 cross p2 in z direction
    
    % Total angular momentum
    L = L1 + L2;
end

