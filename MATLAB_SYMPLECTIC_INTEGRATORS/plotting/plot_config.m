% Configuration for Plotting

% Set the default line width for all future plots
set(0, 'DefaultLineLineWidth', 2);

% Define the new custom color order with darker green and orange
customColorOrder = [
    0 0 1;         % Blue
    1 0 0;         % Red
    0 0.5 0;       % Dark Green
    0.6 0.4 0.2;   % Brown
    1 0.5 0;       % Orange
    1 0.75 0.8;    % Pink
    0 0 0;         % Black
    0.5 0.5 0.5;   % Gray
    0.5 0 0.5;     % Violet
];

% Set the default color order for all future plots
set(0, 'DefaultAxesColorOrder', customColorOrder);

% Set default font properties for all text elements
set(0, 'DefaultAxesFontName', 'Arial');             % Set font to Arial for axes
set(0, 'DefaultAxesFontSize', 12);                  % Set font size to 12 for axes
set(0, 'DefaultAxesFontWeight', 'bold');            % Set font weight to bold for axes
set(0, 'DefaultTextFontName', 'Arial');             % Set font to Arial for text
set(0, 'DefaultTextFontSize', 12);                  % Set font size to 12 for text
set(0, 'DefaultTextFontWeight', 'bold');            % Set font weight to bold for text
set(0, 'DefaultLegendFontName', 'Arial');           % Set font to Arial for legend
set(0, 'DefaultLegendFontSize', 12);                % Set font size to 12 for legend
set(0, 'DefaultLegendFontWeight', 'bold');          % Set font weight to bold for legend

