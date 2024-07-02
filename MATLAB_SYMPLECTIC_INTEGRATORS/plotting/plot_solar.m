function ax = plot_solar(T, Y, folder, name, varargin)
    % Function to plot trajectories of solar system bodies in 3D
    % Inputs:
    %   T: Time vector
    %   Y: Position data matrix (size: [n x l])
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
    if nargin > 5
        video = varargin{1};
        step = varargin{2};
        p = varargin{3};
        ax = varargin{4};
    end

    % Plot Total Trajectories
    fig = figure; 
    legendentries = cell(1, l / 6);
    legendentries{1} = 'Sun';
    legendentries{2} = 'Jupiter';
    legendentries{3} = 'Saturn';
    legendentries{4} = 'Uranus';
    legendentries{5} = 'Neptune';
    legendentries{6} = 'Pluto';

    for i = 1:3:l / 2
        plot3(Y(:, i), Y(:, i + 1), Y(:, i + 2))
        xlabel('x');
        ylabel('y');
        zlabel('z');
        hold on
    end
    axis equal
    legend(legendentries, 'FontSize', 12)
    title('Total trajectory')
    if ax == 0
        ax = axis; % Extract the axis for the video
    else
        axis(ax)
    end
    hold off

    % Specify the folder and filename for saving
    relativeFolderPath = sprintf('results/%s', folder); % Relative folder path (subfolder)
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
        legendentries2 = cell(1, l / 3);
        for i = 1:l / 3
            if mod(i, 2) == 1
                legendentries2{i} = '';
            end
            if mod(i, 2) == 0
                legendentries2{i} = legendentries{i / 2};
            end
        end

        % Define the name and format of the video file
        videoFileName = sprintf('video_%s.mp4', name);
        videoPath = fullfile(relativeFolderPath, videoFileName);
        video = VideoWriter(videoPath, 'MPEG-4');
        video.FrameRate = 10;
        open(video);
        
        fig2 = figure(10);
        for k = 1:step:n
            figure(10)
            km = max(1, k - floor(n / 5));
            titol = strcat('t = ', num2str(T(k)));
            for i = 1:3:l / 2
                currentAxes = gca;
                colorOrder = currentAxes.ColorOrder;
                nextColor = colorOrder(mod(size(currentAxes.Children, 1), size(colorOrder, 1)) + 1, :);
                plot3(Y(k, i), Y(k, i + 1), Y(k, i + 2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', nextColor);
                xlabel('x');
                ylabel('y');
                zlabel('z');
                hold on
                plot3(Y(km:k, i), Y(km:k, i + 1), Y(km:k, i + 2), 'Color', nextColor);
                hold on
            end
            legend(legendentries2, 'FontSize', 12)
            title(titol)
            hold off
            axis(ax)
            axis square
            grid on
            % Capture the plot as a frame
            pause(p)
            frame = getframe(fig2);
            % Write the frame to the video file
            writeVideo(video, frame);
        end
    end
    
end

