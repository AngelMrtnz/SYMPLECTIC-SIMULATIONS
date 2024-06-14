function n_body_2d
	clc
	clear all
	close all
    addpath('examples');
    addpath('methods');
    plot_config;

	m=[1/3,1/3,1/3];
    G=1;

    c1=  0.1546025015193900;
    c2= -0.0987561640960160;
    c3= -1.1843704912618038;
    c4= 0.2521819944674475;

    %c1=  0.1540289328762600;
    %c2=-0.0932474525425240;
    %c3=0.9635026027520550;
    %c4=0.5076890556787989;

    ob10=[-2*c1,0];
    ob20=[c1,c2];
    ob30=[c1,-c2];
    p_ob10=[m(1)*0,m(1)*-2*c4];
    p_ob20=[m(2)*c3,m(2)*c4];
    p_ob30=[m(3)*-c3,m(3)*c4];

	ci = [ob10, ob20,ob30, p_ob10, p_ob20,p_ob30];
	ti = 0;
	tf = 2000;
	h = .001;

    % Start the timer for function 1
    tic;
    [T1,Y1] = method_sv_faster(@(t,y) f1(t,y,m,G),@(t,y) f2(t,y,m,G), [ti:h:tf], ci);
    elapsedTime1 = toc;
    disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']); 
    H1=Hamiltonian(T1,Y1,m,G);
    H1_rel=H1-H1(1);
    L1=angular_momentum(T1,Y1,m,G);
    options = odeset('RelTol',1e-12);
    % Start the timer for function 2
    
    tic;
    [T2,Y2] = method_rk3(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime2 = toc;
    disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']); 
    H2=Hamiltonian(T2,Y2,m,G);
    H2_rel=H2-H2(1);
    L2=angular_momentum(T2,Y2,m,G);
    
    
    % Start the timer for function 3
    tic;
    [T3,Y3] = method_rk5(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime3 = toc;
    disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']); 
    H3=Hamiltonian(T3,Y3,m,G);
    H3_rel=H3-H3(1);
    L3=angular_momentum(T3,Y3,m,G);
    
    % Start the timer for function 4
    tic;
    [T4,Y4] = method_midpoint(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime4 = toc;
    disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']); 
    H4=Hamiltonian(T4,Y4,m,G);
    H4_rel=H4-H4(1);
    L4=angular_momentum(T4,Y4,m,G);

    % Start the timer for function 5
    tic;
    [T5,Y5] = method_ralston(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
    elapsedTime5 = toc;
    disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']); 
    H5=Hamiltonian(T5,Y5,m,G);
    H5_rel=H5-H5(1);
    L5=angular_momentum(T5,Y5,m,G);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% ENERGY AND ANGULAR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig=figure(3);
    set(fig, 'Position', [100, 100, 1000, 500]); % [left, bottom, width, height]

    subplot(1,2,1)
    plot(T1,H1_rel,T2,H2_rel,T3,H3_rel,":", T4,H4_rel,"-.",T5,H5_rel,'--')
    title('Difference in Energy')
    xlabel('time');
    ylabel('H-H(0)');
    %legend('S-V 1','S-V 2','R-K 5','FontSize', 10,'Location', 'northwest')
    legend('S-V 2','RK-3','RK-5','Midpoint','Ralston','FontSize', 10,'Location', 'northwest')
    subplot(1,2,2)
    plot(T1, (L1), T2, (L2),'--',T3, (L3),':', T4,(L4),'-.',T5, (L5))
    title('Difference angular momentum')
    xlabel('time');
    ylabel('L-L(0) (norm)');
    %legend('S-V 1','S-V 2','R-K 5','FontSize', 10,'Location', 'northwest')
    legend('S-V 2','RK-3','RK-5','Midpoint','Ralston','FontSize', 10,'Location', 'northwest')
    %legend('S-V 1','Heun','Midpoint','Ralston','RK-5','FontSize', 10,'Location', 'northwest')

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
    relativeFolderPath = 'results/three_body/test'; % Relative folder path (subfolder)
    fileName = 'energy_3_8.pdf';
    
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
    ax = plot_n_bodies_2d(T5,Y5,'three_body/test','3_ralston8');
    plot_n_bodies_2d(T4,Y4,'three_body/test','3_midpoint8',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T3,Y3,'three_body/test','3_rk38',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T2,Y2,'three_body/test','3_rk58',0,floor((tf/h)/300),.001,ax);
    plot_n_bodies_2d(T1,Y1,'three_body/test','3_sv218',0,floor((tf/h)/300),.001,ax);
    

end
function dy=f1(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/4;
    dy=zeros(l/2,1);
    aux=1; %we only have N masses and so we need to add 1 for every three indices
    for i=1:l/2
        dy(i)=y(i+l/2)/m(aux);
        if  mod(i,2)==0 
            aux=aux+1;
        end
    end
end
function dy=f2(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/4;
    dy=zeros(l/2,1);
    x_c=[]; %Coordinate vectors, i.e x_c=(x1,x2,x3,...,x_n)
    y_c=[];
    r=zeros(N,N); %matrix with all r_ij, relative distances
    for i=1:N
        ob1=y(2*i-1:2*i);
        for j=1:i-1
            ob2=y(2*j-1:2*j);

            r(i,j) = dist(ob1,ob2);
        end
    end
    r=r+r'; %it is symmetric

    aux=1; %we only have N masses and so we need to add 1 for every two indices
    for i=1:l/2
      %filling in the coordinate vectors and updating j;
        if  mod(i,2)==1 
            x_c=[x_c y(i)];
        end
        if  mod(i,2)==0 
            y_c=[y_c y(i)];
        end
    end
    
    aux=0;
    for i=l/2+1:2:l %loop for the derivate of moments
        j=i-(2*N)-(aux);
         for k=1:N
           if k~=j
               dy(i-l/2) = dy(i-l/2) -G*m(j)*m(k)*(x_c(j)-x_c(k))/r(j,k)^3;
               dy(i-l/2+1) = dy(i-l/2+1) -G*m(j)*m(k)*(y_c(j)-y_c(k))/r(j,k)^3;
           end  
         end
         aux=aux+1;
    end
    
end


function dy=f(t,y,m,G)
    sup=size(y);
    l=sup(1)*sup(2);
    N=l/4; 
    dy=zeros(l,1);
    x_c=[]; %Coordinate vectors, i.e x_c=(x1,x2,...,x_N)
    y_c=[];
    r=zeros(N,N); %matrix with all r_ij, relative distances
    for i=1:N
        ob1=y(2*i-1:2*i);
        for j=1:i-1
            ob2=y(2*j-1:2*j);

            r(i,j) = dist(ob1,ob2);
        end
    end
    r=r+r'; %it is symmetric

    aux=1; %we only have N masses and so we need to add 1 for every two indices
    for i=1:l/2
        dy(i)=y(i+l/2)/m(aux);
        %filling in the coordinate vectors and updating j;
        if  mod(i,2)==1 
            x_c=[x_c y(i)];
        end
        if  mod(i,2)==0 
            y_c=[y_c y(i)];
            aux=aux+1;
        end
    end
    
    aux=0;
    for i=l/2+1:2:l %loop for the derivate of moments
        j=i-(2*N)-(aux);
         for k=1:N
           if k~=j
               dy(i) = dy(i) -G*m(j)*m(k)*(x_c(j)-x_c(k))/r(j,k)^3;
               dy(i+1) = dy(i+1) -G*m(j)*m(k)*(y_c(j)-y_c(k))/r(j,k)^3;
           end  
         end
         aux=aux+1;
         
    end
    
end

function H=Hamiltonian(t,y,m,G)
    N=length(m);
    l=N*4;
    H=0;
    for i=1:N
        for j=1:i-1
            ob1=[y(:,2*j-1:2*j)];
            ob2=[y(:,2*i-1:2*i)];
            H=H-G*m(i)*m(j)*(dist2(ob1,ob2).^(-1));
        end

    end
    
        aux=1; %we only have N masses and so we need to add 1 for every two indices
    for i=1:l/2
       H=H+y(:,i+l/2).^2/(2*m(aux));
       if  mod(i,2)==0 
            aux=aux+1;
        end
    end

end

function L_mag = angular_momentum(t, y, m,G)
    % Parameters:
    % t - time vector
    % y - state matrix [positions (2N) and velocities (2N) at each time step]
    % m - masses of the bodies (Nx1 vector)
    
    N = length(m);  % Number of bodies
    l = N * 4;  % Length of state vector for one time step
    num_steps = length(t);  % Number of time steps
    
    % Initialize angular momentum magnitude vector
    L_mag = zeros(num_steps, 1);
    
    % Loop over each time step to compute the angular momentum
    for k = 1:num_steps
        % Initialize angular momentum vector for the current time step
        L = zeros(1, 2);
        
        % Loop over each body to compute the angular momentum
        for i = 1:N
            % Extract position and velocity for the current body
            r = y(k, 2*i-1:2*i);  % Position vector
            p = y(k, l/2 + 2*i-1:l/2 + 2*i);  % Momentum vector
            r_extended=[r,0]';
            p_extended=[p,0]';
            
            % Compute the angular momentum for the current body and add to total
            L = L + cross(r_extended,p_extended);
        end
        
        % Compute the magnitude of the angular momentum vector for the current time step
        L_mag(k) = abs(L(3)-L_mag(1));
    end
end


function rij=dist(ob1,ob2)
	[fil,col] = size(ob1);
	if (col < fil)
		ob1 = ob1';
        ob2=ob2';
	end
    rij=((ob1(:,1)-ob2(:,1)).^2+(ob1(:,2)-ob2(:,2)).^2).^(1/2);
end

function rij=dist2(ob1,ob2)
    rij=((ob1(:,1)-ob2(:,1)).^2+(ob1(:,2)-ob2(:,2)).^2).^(1/2);
end
