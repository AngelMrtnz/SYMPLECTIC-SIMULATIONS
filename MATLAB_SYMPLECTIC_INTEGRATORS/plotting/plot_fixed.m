function ax = plot_fixed(T, Y, folder, name, varargin)
    % Function to plot fixed points of trajectories
    % Inputs:
    %   T: Time vector
    %   Y: Position and momenta data matrix 
    %   folder: Folder name for saving results
    %   name: Name identifier for the plot
    %   varargin: Optional arguments for video generation
    % Outputs:
    %   ax: Axis limits used for the plots

    [n, l] = size(Y); % Get the dimensions of Y

    % Initialize parameters for plotting
    video = 0; % Flag for generating video (0: no, 1: yes)
    p = 0.001; % Pause time between frames
    step = floor(1 / n); % Step size for frames
    ax = 0; % Axis limits

    % Check if additional arguments are provided
    if nargin == 5
        video = varargin{1};
    end

    if nargin >= 5
        video = varargin{1};
        step = varargin{2};
        p = varargin{3};
        ax = varargin{4};
    end

    % Plot Total Trajectories
    fig = figure; 
    for i = 1:2:l / 2 - 2
        currentAxes = gca;
        colorOrder = currentAxes.ColorOrder;
        nextColor = colorOrder(mod(size(currentAxes.Children, 1), size(colorOrder, 1)) + 1, :);
        plot(Y(1, i), Y(1, i + 1), 'o', 'MarkerSize', 5, 'MarkerFaceColor', nextColor)
        hold on
    end
    plot(Y(:, l / 2 - 1), Y(:, l / 2))
    xlabel('x');
    ylabel('y');
    title('Total trajectory')
    axis equal
    if ax == 0
        ax = axis; % Extract the axis for the video
    else
        axis(ax)
    end
    hold off

    % Specify the folder and filename for saving
    relativeFolderPath =  sprintf('results/%s', folder);
    fileName = sprintf('trajectory_%s.pdf', name);

    % Ensure the folder exists, if not, create it
    if ~exist(relativeFolderPath, 'dir')
        mkdir(relativeFolderPath);
    end
    
    % Full file path
    filePath = fullfile(relativeFolderPath, fileName);
    
    % Adjust the size of the paper to match the figure
    set(fig, 'Units', 'Inches');
    pos = get(fig, 'Position');
    set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    
    % Save the figure as a pdf
    saveas(fig, filePath);

    % Plot Movement (Generate Video)
    if video == 1
        videoFileName = sprintf('video_%s.mp4', name);
        videoPath = fullfile(relativeFolderPath, videoFileName);
        video = VideoWriter(videoPath, 'MPEG-4');
        video.FrameRate = 10;
        open(video);
    
        fig2 = figure(2);
        for k = 1:step:n
            figure(2)
            km = max(1, k - step * 2);
            titol = strcat('t = ', num2str(T(k)));
            for i = 1:2:(l / 2)-2
                plot(Y(1, i), Y(1, i + 1), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
                hold on;
            end
            plot(Y(km:k, l/2-1), Y(km:k, l/2),'LineWidth',5);
            xlabel('x');
            ylabel('y');
            title(titol)
            hold off
            axis(ax)
            axis square
            grid on
            pause(p)
            frame = getframe(fig2);
            writeVideo(video, frame);
        end
    end
end


