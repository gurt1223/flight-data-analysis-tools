function formatPlot()

% Format Plot for Consistent Visualization
%
% DESCRIPTION:
%   This function sets default plotting parameters such as font size,
%   grid lines, tick marks, and line width to ensure consistent visualization
%   of plots.
%
% INPUTS:
%   None
%
% OUTPUTS:
%   None (Modifies default plotting settings)
% 
% WRITTEN BY:
%   Garrett D. Asper
%   Virginia Tech
%   Email: garrettasper@vt.edu
%
% HISTORY:
%   03 JUL 2024 - Created and debugged, GDA
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',14);
set(0,'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on')   % Turn grid on for all plots
set(0, 'DefaultAxesYMinorTick', 'on');                      % Turn on minor tick marks for y-axis
set(0, 'DefaultAxesYTickMode', 'auto');                     % or 'manual' if you want to specify ticks manually
set(0,'DefaultAxesFontSize', 15)                            % Change the default figure axes font size
set(0, 'DefaultLineLineWidth', 1.1)                         % Set default plot line width
set(0, 'DefaultAxesLineWidth', 1.1)                         % Set default axes line width

end