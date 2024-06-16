% This script simulates the motion of fixed n-body systems using different numerical methods and compares their performance in terms of computational time and energy conservation.
function n_body_2d_fixed
% Clear workspace
clc
clear all
close all

% Add paths to the required folders
addpath('methods');
addpath('plotting')
addpath('examples');

% Load plotting configuration
plot_config;

% Define the masses of the bodies
m = [100,100,100,1];

% Define gravitational constant
G = 1;

% Define initial positions of the bodies
ob10 = [5,0];
ob20 = [0,0];
ob30 = [10,0];
ob40 = [-5,0];
p_ob10 = [0,0];
p_ob20 = [0,0];
p_ob30 = [0,0];
p_ob40 = [0,-7];

% Combine initial positions into a single array
ci = [ob10, ob20,ob30,ob40, p_ob10, p_ob20,p_ob30,p_ob40];

% Define simulation parameters
ti = 0;             % Initial time
tf = 300;          % Final time
h = .005;           % Step size
options = odeset('RelTol',1e-12);  % ODE options

% Perform simulations using different numerical methods and measure elapsed time
tic;
[T1,Y1] = method_sv_faster2(@(t,y) f1(t,y,m,G),@(t,y) f2(t,y,m,G), [ti:h:tf], ci);
elapsedTime1 = toc;
disp(['Elapsed time for function 1: ', num2str(elapsedTime1), ' seconds']);

tic;
[T2,Y2] = method_rk3(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
elapsedTime2 = toc;
disp(['Elapsed time for function 2: ', num2str(elapsedTime2), ' seconds']); 

tic;
[T3,Y3] = method_rk5(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
elapsedTime3 = toc;
disp(['Elapsed time for function 3: ', num2str(elapsedTime3), ' seconds']); 

tic;
[T4,Y4] = method_midpoint(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
elapsedTime4 = toc;
disp(['Elapsed time for function 4: ', num2str(elapsedTime4), ' seconds']); 

tic;
[T5,Y5] = method_ralston(@(t,y) f(t,y,m,G), [ti:h:tf], ci);
elapsedTime5 = toc;
disp(['Elapsed time for function 5: ', num2str(elapsedTime5), ' seconds']); 

% Calculate relative Hamiltonians for each method
H1 = Hamiltonian(T1,Y1,m,G);
H1_rel = H1 - H1(1);

H2 = Hamiltonian(T2,Y2,m,G);
H2_rel = H2 - H2(1);

H3 = Hamiltonian(T3,Y3,m,G);
H3_rel = H3 - H3(1);

H4 = Hamiltonian(T4,Y4,m,G);
H4_rel = H4 - H4(1);

H5 = Hamiltonian(T5,Y5,m,G);
H5_rel = H5 - H5(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ENERGY AND ANGULAR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot relative energy over time for each method
fig = figure(3);
set(fig, 'Position', [100, 100, 1000, 500]); % [left, bottom, width, height]
plot(T1,H1_rel,T2,H2_rel,'--',T3,H3_rel,":", T4,H4_rel,"-.",T5,H5_rel)
title('Energy (relative)')
xlabel('time');
ylabel('H-H(0)');
legend('S-V 1','RK-3','RK-5','Midpoint','Ralston','FontSize', 10,'Location', 'northwest')

% Specify the relative folder and filename for saving the figure
relativeFolderPath = 'results/n_fixed_body'; % Relative folder path (subfolder)
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

% Plot the fixed n-body systems for each method
    ax=plot_fixed(T3,Y3,'n_fixed_body','2_rk5');
    plot_fixed(T4,Y4,'n_fixed_body','2_midpoint',0,floor((tf/h)/300),.001,ax);
    plot_fixed(T5,Y5,'n_fixed_body','2_ralston',0,floor((tf/h)/300),.001,ax);
    plot_fixed(T2,Y2,'n_fixed_body','2_rk3',0,floor((tf/h)/300),.001,ax);
    plot_fixed(T1,Y1,'n_fixed_body','2_sv21',1,10,.001,ax);

end

% Function to compute derivatives for the first function
function dy = f1(t, y, m, G)
    sup = size(y);
    l = sup(1) * sup(2);
    N = l / 4;
    dy = zeros(l / 2, 1);
    % Fix n-1 bodies
    dy(l / 2 - 1) = y(l - 1) / m(N);
    dy(l / 2) = y(l) / m(N);
end

% Function to compute derivatives for the second function
function dy = f2(t, y, m, G)
    sup = size(y);
    l = sup(1) * sup(2);
    N = l / 4;
    dy = zeros(l / 2, 1);
    x_c = []; % Coordinate vectors, i.e x_c=(x1,x2,...,x_N)
    y_c = [];
    r = zeros(N, N); % Matrix with all r_ij, relative distances
    i = N;
    ob1 = y(2 * i - 1:2 * i);
    for j = 1:i - 1
        ob2 = y(2 * j - 1:2 * j);
        r(i, j) = dist(ob1, ob2);
    end
    r = r + r'; % It is symmetric
    for i = 1:l / 2
        % Filling in the coordinate vectors and updating j;
        if mod(i, 2) == 1
            x_c = [x_c y(i)];
        end
        if mod(i, 2) == 0
            y_c = [y_c y(i)];
        end
    end
    % Fix n-1 bodies
    i = l - 1;
    j = i - (2 * N) - (N - 1);
    for k = 1:N
        if k ~= j
            dy(i - l / 2) = dy(i - l / 2) - G * m(j) * m(k) * (x_c(j) - x_c(k)) / r(j, k)^3;
            dy(i + 1 - l / 2) = dy(i + 1 - l / 2) - G * m(j) * m(k) * (y_c(j) - y_c(k)) / r(j, k)^3;
        end
    end
end

% Function to compute derivatives for the third function
function dy = f(t, y, m, G)
    sup = size(y);
    l = sup(1) * sup(2);
    N = l / 4;
    dy = zeros(l, 1);
    x_c = []; % Coordinate vectors, i.e x_c=(x1,x2,...,x_N)
    y_c = [];
    r = zeros(N, N); % Matrix with all r_ij, relative distances
    i = N;
    ob1 = y(2 * i - 1:2 * i);
    for j = 1:i - 1
        ob2 = y(2 * j - 1:2 * j);
        r(i, j) = dist(ob1, ob2);
    end
    r = r + r'; % It is symmetric
    for i = 1:l / 2
        % Filling in the coordinate vectors and updating j;
        if mod(i, 2) == 1
            x_c = [x_c y(i)];
        end
        if mod(i, 2) == 0
            y_c = [y_c y(i)];
        end
    end
    % Fix n-1 bodies
    dy(l / 2 - 1) = y(l - 1) / m(N);
    dy(l / 2) = y(l) / m(N);
    i = l - 1;
    j = i - (2 * N) - (N - 1);
    for k = 1:N
        if k ~= j
            dy(i) = dy(i) - G * m(j) * m(k) * (x_c(j) - x_c(k)) / r(j, k)^3;
            dy(i + 1) = dy(i + 1) - G * m(j) * m(k) * (y_c(j) - y_c(k)) / r(j, k)^3;
        end
    end
end

% Function to compute Hamiltonian
function H = Hamiltonian(t, y, m, G)
    N = length(m);
    l = N * 4;
    H = 0;
    i = N;
    ob2 = [y(:, 2 * i - 1:2 * i)];
    for j = 1:i - 1
        ob1 = [y(:, 2 * j - 1:2 * j)];
        H = H - G * m(i) * m(j) * (dist2(ob1, ob2).^(-1));
    end
    for i = l / 2 - 1:l / 2
        H = H + y(:, i + l / 2).^2 / (2 * m(N));
    end
end

% Function to compute distance between two points
function rij = dist(ob1, ob2)
    [fil, col] = size(ob1);
    if (col < fil)
        ob1 = ob1';
        ob2 = ob2';
    end
    rij = ((ob1(:, 1) - ob2(:, 1)).^2 + (ob1(:, 2) - ob2(:, 2)).^2).^(1 / 2);
end

% Function to compute square of distance between two points
function rij = dist2(ob1, ob2)
    rij = ((ob1(:, 1) - ob2(:, 1)).^2 + (ob1(:, 2) - ob2(:, 2)).^2).^(1 / 2);
end

