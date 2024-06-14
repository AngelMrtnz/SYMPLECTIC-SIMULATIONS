% Simulate an n-body system with different numerical integration methods
function n_body
    clc; % Clear command window
    clear all; % Clear workspace
    close all; % Close all figures

    % Add necessary paths
    addpath('methods');
    addpath('examples');
    addpath('plotting')

    % Configure plot settings
    plot_config;
    
    % Load parameters for the n-body example
    [G,m,ci] = n_body_examples('case 1');
    
    % Define simulation time parameters
    ti = 0;
    tf = 50;
    h = 0.001;

    % Start the timer for function 1
    tic;
    [T1,Y1] = method_sv_faster2(@(t,y) f1(t,y,m,G),@(t,y) f2(t,y,m,G), [ti:h:tf], ci);
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']); 
    H1 = Hamiltonian(T1,Y1,m,G);
    H1_rel = H1 - H1(1);
    L1 = angular_momentum(T1,Y1,m,G);
    
    % Start the timer for function 2
    tic;
    [T2,Y2] = method_rk3(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']); 
    H2 = Hamiltonian(T2,Y2,m,G);
    H2_rel = H2 - H2(1);
    L2 = angular_momentum(T2,Y2,m,G);
    
    % Start the timer for function 3
    tic;
    [T3,Y3] = method_rk5(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']); 
    H3 = Hamiltonian(T3,Y3,m,G);
    H3_rel = H3 - H3(1);
    L3 = angular_momentum(T3,Y3,m,G);
    
    % Start the timer for function 4
    tic;
    [T4,Y4] = method_midpoint(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']); 
    H4 = Hamiltonian(T4,Y4,m,G);
    H4_rel = H4 - H4(1);
    L4 = angular_momentum(T4,Y4,m,G);

    % Start the timer for function 5
    tic;
    [T5,Y5] = method_ralston(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']); 
    H5 = Hamiltonian(T5,Y5,m,G);
    H5_rel = H5 - H5(1);
    L5 = angular_momentum(T5,Y5,m,G);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% ENERGY AND ANGULAR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create a new figure for plotting
    fig = figure(3);
    set(fig, 'Position', [100, 100, 1000, 500]); % [left, bottom, width, height]

    % Plot difference in energy
    subplot(1,2,1)
    plot(T1, H1_rel, T2, H2_rel, '--', T3, H3_rel, ':', T4, H4_rel, '-.', T5, H5_rel)
    title('Difference in energy')
    xlabel('time');
    ylabel('H-H(0)');
    legend('S-V 2','RK-3','RK-5','Midpoint','Ralston','FontSize', 10,'Location', 'northwest')
    
    % Plot difference in angular momentum
    subplot(1,2,2)
    plot(T1, L1, T2, L2, '--', T3, L3, ':', T4, L4, '-.', T5, L5)
    title('Difference angular momentum')
    xlabel('time');
    ylabel('L-L(0) (norm)');
    legend('S-V 2','RK-3','RK-5','Midpoint','Ralston','FontSize', 10,'Location', 'northwest')

    % Adjust the position of the subplots
    h1 = subplot(1, 2, 1);
    pos1 = get(h1, 'Position');
    pos1(1) = pos1(1) - 0.06;
    pos1(3) = pos1(3) * 1.1;
    set(h1, 'Position', pos1);
    
    h2 = subplot(1, 2, 2);
    pos2 = get(h2, 'Position');
    pos2(1) = pos2(1) + 0;
    pos2(3) = pos2(3) * 1.1;
    set(h2, 'Position', pos2);

    % Specify the relative folder and filename for saving results
    relativeFolderPath = 'results/n_body';
    fileName = 'energy_3_test.pdf';
    
    % Ensure the folder exists, if not, create it
    if ~exist(relativeFolderPath, 'dir')
        mkdir(relativeFolderPath);
    end
    
    % Full file path for saving the figure
    filePath = fullfile(relativeFolderPath, fileName);
    
    % Adjust the size of the paper to match the figure
    set(fig, 'Units', 'Inches');
    pos = get(fig, 'Position');
    set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Save the figure as a PDF file
    saveas(fig, filePath);
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot 3D trajectories for each integration method
    ax = plot_n_bodies_3d(T3,Y3,'n_body','3_rk5');
    plot_n_bodies_3d(T4,Y4,'n_body','3_midpoint',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_3d(T5,Y5,'n_body','3_ralston',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_3d(T2,Y2,'n_body','3_rk3',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_3d(T1,Y1,'n_body','3_sv21',1,floor((tf/h)/300),.001,ax);
end

%Differential equations
function dy=f1(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/6;
    dy=zeros(l/2,1);
    aux=1; %we only have N masses and so we need to add 1 for every three indices
    for i=1:l/2
        dy(i)=y(i+l/2)/m(aux);
        if  mod(i,3)==0 
            aux=aux+1;
        end
    end
end

function dy=f2(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/6;
    dy=zeros(l/2,1);
    x_c=[]; %Coordinate vectors, i.e x_c=(x1,x2,x3,...,x_n)
    y_c=[];
    z_c=[];
    r=zeros(N,N); %matrix with all r_ij, relative distances
    for i=1:N
        for j=1:i-1
            ob1=y(3*i-2:3*i);
            ob2=y(3*j-2:3*j);

            r(i,j) = dist(ob1,ob2);
        end
    end
    r=r+r'; %it is symmetric

    aux=1; %we only have N masses and so we need to add 1 for every three indices
    for i=1:l/2
      %filling in the coordinate vectors and updating j;
        if  mod(i,3)==1 
            x_c=[x_c y(i)];
        end
        if  mod(i,3)==2 
            y_c=[y_c y(i)];
        end
        if  mod(i,3)==0 
            z_c=[z_c y(i)];
            aux=aux+1;
        end
    end
    
    aux=0;
    for i=l/2+1:3:l %loop for the derivate of moments
        j=i-(3*N)-2*(aux);
         for k=1:N
           if k~=j
               dy(i-l/2) = dy(i-l/2) -G*m(j)*m(k)*(x_c(j)-x_c(k))/r(j,k)^3;
               dy(i-l/2+1) = dy(i-l/2+1) -G*m(j)*m(k)*(y_c(j)-y_c(k))/r(j,k)^3;
               dy(i-l/2+2) = dy(i-l/2+2) -G*m(j)*m(k)*(z_c(j)-z_c(k))/r(j,k)^3;
           end  
         end
         aux=aux+1;
    end
    
end

function dy=f(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/6;
    dy=zeros(l,1);
    x_c=[]; %Coordinate vectors, i.e x_c=(x1,x2,x3,...,x_n)
    y_c=[];
    z_c=[];
    r=zeros(N,N); %matrix with all r_ij, relative distances
    for i=1:N
        for j=1:i-1
            ob1=y(3*i-2:3*i);
            ob2=y(3*j-2:3*j);

            r(i,j) = dist(ob1,ob2);
        end
    end
    r=r+r'; %it is symmetric

    aux=1; %we only have N masses and so we need to add 1 for every three indices
    for i=1:l/2
        dy(i)=y(i+l/2)/m(aux);
        %filling in the coordinate vectors and updating j;
        if  mod(i,3)==1 
            x_c=[x_c y(i)];
        end
        if  mod(i,3)==2 
            y_c=[y_c y(i)];
        end
        if  mod(i,3)==0 
            z_c=[z_c y(i)];
            aux=aux+1;
        end
    end
    

    aux=0;
    for i=l/2+1:3:l %loop for the derivate of moments
        j=i-(3*N)-2*(aux);
         for k=1:N
           if k~=j
               dy(i) = dy(i) -G*m(j)*m(k)*(x_c(j)-x_c(k))/r(j,k)^3;
               dy(i+1) = dy(i+1) -G*m(j)*m(k)*(y_c(j)-y_c(k))/r(j,k)^3;
               dy(i+2) = dy(i+2) -G*m(j)*m(k)*(z_c(j)-z_c(k))/r(j,k)^3;
           end  
         end
         aux=aux+1;
    end
    
end

%Hamiltonian function
function H=Hamiltonian(t,y,m,G)
    N=length(m);
    l=N*6;
    H=0;
    for i=1:N
        for j=1:i-1
            ob1=[y(:,3*j-2:3*j)];
            ob2=[y(:,3*i-2:3*i)];
            H=H-G*m(i)*m(j)*(dist2(ob1,ob2).^(-1));
        end

    end
    
        aux=1; %we only have N masses and so we need to add 1 for every three indices
    for i=1:l/2
       H=H+y(:,i+l/2).^2/(2*m(aux));
       if  mod(i,3)==0 
            aux=aux+1;
        end
    end

end
%norm of Angular momentum
function L_mag = angular_momentum(t, y, m,G)
    % Parameters:
    % t - time vector
    % y - state matrix [positions (3N) and velocities (3N) at each time step]
    % m - masses of the bodies (Nx1 vector)
    
    N = length(m);  % Number of bodies
    l = N * 6;  % Length of state vector for one time step
    num_steps = length(t);  % Number of time steps
    
    % Initialize angular momentum magnitude vector
    L_mag = zeros(num_steps, 1);
    
    % Loop over each time step to compute the angular momentum
            % Initialize angular momentum vector for the current time step
        L_ini = zeros(1, 3);
        
        % Loop over each body to compute the angular momentum
        for i = 1:N
            % Extract position and velocity for the current body
            r = y(1, 3*i-2:3*i);  % Position vector
            p = y(1, l/2 + 3*i-2:l/2 + 3*i);  % Momentum vector
            
            
            % Compute the angular momentum for the current body and add to total
            L_ini = L_ini + cross(r, p);
        end
    for k = 2:num_steps
        % Initialize angular momentum vector for the current time step
        L = zeros(1, 3);
        
        % Loop over each body to compute the angular momentum
        for i = 1:N
            % Extract position and velocity for the current body
            r = y(k, 3*i-2:3*i);  % Position vector
            p = y(k, l/2 + 3*i-2:l/2 + 3*i);  % Momentum vector
            
            
            % Compute the angular momentum for the current body and add to total
            L = L + cross(r, p);
        end
        
        % Compute the magnitude of the angular momentum vector for the current time step
        L_mag(k) = norm(L-L_ini);
    end
end

%Functions for distance between two points
function rij=dist(ob1,ob2)
	[fil,col] = size(ob1);
	if (col < fil)
		ob1 = ob1';
        ob2=ob2';
	end
    rij=((ob1(:,1)-ob2(:,1)).^2+(ob1(:,2)-ob2(:,2)).^2+(ob1(:,3)-ob2(:,3)).^2).^(1/2);
end

function rij=dist2(ob1,ob2)
    rij=((ob1(:,1)-ob2(:,1)).^2+(ob1(:,2)-ob2(:,2)).^2+(ob1(:,3)-ob2(:,3)).^2).^(1/2);
end

