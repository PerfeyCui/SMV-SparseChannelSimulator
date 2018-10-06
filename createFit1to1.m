function [pd1] = createFit1to1(arg_1,fitType,plotType)
% Input arg_1:input data;
        %fitType: distribution types. 'normal','lognormal' etc
        %PlotType: plot show type. 'cdf','pdf','Prob'(Probability figures);
%CREATEFIT    Create plot of datasets and fits
%   [PD1,PD2] = CREATEFIT(ARG_1)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
  %改了两处，为了同图画拟合曲线；
% Force all inputs to be column vectors
arg_1 = arg_1(:);
% initializations:
if nargin <2, 
    fitType = 'norm';
    plotType = 'cdf';
elseif nargin <3
    plotType = 'cdf';
end

% Prepare figure
%clf; %%%%%%%%%%%%%%%%%%%%%%
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "Tap1Loc(6,:) data"
switch lower(plotType)
    case{'cdf'}
%cdf
[CdfY,CdfX] = ecdf(arg_1,'Function','cdf');  % compute empirical function
hLine = stairs(CdfX,CdfY,'Color',[0.333333 0 0.666667],'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Cumulative probability')
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Tap1Loc(6,:) data';
    case{'pdf'}
%pdf
[CdfF,CdfX] = ecdf(arg_1,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(arg_1,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
    case{'pro','probability','probabilityplot'}
%probbility plot
hLine = probplot('normal',arg_1,[],[],'noref');
set(hLine,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Data');
ylabel('Probability')
    otherwise
    %probbility plot
hLine = probplot('normal',arg_1,[],[],'noref');
set(hLine,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Data');
ylabel('Probability')
end
LegHandles(end+1) = hLine;
LegText{end+1} = [ 'Original data'];    
% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);



% --- Create fit "LogNormFit"
% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('lognormal',[ -1.449972153488, 0.485820282068])
pd1 = fitdist(arg_1,fitType);%'lognormal'
switch lower(plotType)
    case{'cdf'}
YPlot = cdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
    case{'pdf'}
YPlot = pdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
    case{'pro','probability','probabilityplot'}
hLine = probplot(gca,pd1);
set(hLine,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2);
    otherwise
       hLine = probplot(gca,pd1);
set(hLine,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2); 
end
LegHandles(end+1) = hLine;
LegText{end+1} = [fitType 'Fit'];%'LogNormFit';

% Adjust figure
box on;
grid on;
%hold off;%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northwest');
set(hLegend,'Interpreter','none');
