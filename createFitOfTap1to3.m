function [pd1,pd2,pd3] = createFitOfTap1to3(arg_1,tag)
%CREATEFIT    Create plot of datasets and fits
%   [PD1,PD2,PD3] = CREATEFIT(ARG_1)
% tag = 'Noo'; turn off pd2 and pd3 fitness;
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  3
%
%   See also FITDIST.

% This function was automatically generated on 12-Jul-2018 23:53:05

% Output fitted probablility distributions: PD1,PD2,PD3

% Data from dataset "TapLevInv(12,:) data":
%    Y = arg_1 (originally TapLevInv(12,:))

% Force all inputs to be column vectors
arg_1 = arg_1(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};

if nargin<2
    tag = 'Yes';
end

% --- Plot data originally in dataset "TapLevInv(12,:) data"
[CdfY,CdfX] = ecdf(arg_1,'Function','cdf');  % compute empirical function
hLine = stairs(CdfX,CdfY,'Color',[0.333333 0 0.666667],'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Cumulative probability')
LegHandles(end+1) = hLine;
LegText{end+1} = 'TapLevInv(12,:) data';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% --- Create fit "HalfNormalFit"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('half normal',[ 1, 2.392697222801])
pd1 = fitdist(arg_1, 'half normal');
YPlot = cdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'HalfNormalFit';

% --- Create fit "InverseGaussian"
if tag =='Yes'
% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd2 = ProbDistUnivParam('inverse gaussian',[ 2.765, 5.958184962317])
pd2 = fitdist(arg_1, 'inverse gaussian');
YPlot = cdf(pd2,XGrid);
hLine = plot(XGrid,YPlot,'Color',[0 0 1],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'InverseGaussian';

% --- Create fit "Poisson"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd3 = ProbDistUnivParam('poisson',[ 2.765])
pd3 = fitdist(arg_1, 'poisson');
YPlot = cdf(pd3,XGrid);
hLine = plot(XGrid,YPlot,'Color',[0.666667 0.333333 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'Poisson';
else
    pd2 = [];pd3=[];
end

% Adjust figure
box on;
grid on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northwest');
set(hLegend,'Interpreter','none');
