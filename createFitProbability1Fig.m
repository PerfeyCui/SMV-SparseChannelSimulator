function [PD] = createFitProbability1Fig(Arg,FitType,OuterFig,col) %arg_1,arg_2
%CREATEFIT    Create plot of datasets and fits
%   [PD1,PD2] = CREATEFIT(ARG_1,ARG_2),
% OuterFig = 'Yes', 利用外在的plotfig；,Noo self define；
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%  %
%   Number of datasets:  2
%   Number of fits:  2
%
%   See also FITDIST.

% This function was automatically generated on 13-Jul-2018 22:45:45

% Output fitted probablility distributions: PD1,PD2

% Data from dataset "TapLoc(5,:) data":
%    Y = arg_1 (originally TapLoc(5,:))

% Data from dataset "TapLoc(12,:) data":
%    Y = arg_2 (originally TapLoc(12,:))

% % Force all inputs to be column vectors
% arg_1 = arg_1(:);
% arg_2 = arg_2(:);

% Prepare figure
%clf;
if nargin <4, col = 'r'; end

if nargin<2,
    FitType = 'lognormal';
    OuterFig = 'Noo';
elseif nargin<3,
    OuterFig = 'Noo';
end

hold on;
if OuterFig == 'Noo'
probplot('normal'); % create empty plot of desired type
title('');
end

if size(Arg,1)==1
   Arg = Arg'; 
end

PD = [];
for ii = 1:size(Arg,2)
    arg_1 = Arg(:,ii);
% --- Plot data originally in dataset "TapLoc(5,:) data"
hLine = probplot(gca,arg_1,[],[],'noref'); % add data to existing plot
set(hLine, 'MarkerSize',6,'Color',col);
xlabel('Data');
ylabel('Probability')

% % --- Plot data originally in dataset "TapLoc(12,:) data"
% hLine = probplot(gca,arg_2,[],[],'noref'); % add data to existing plot
% set(hLine,'Color',[0.333333 0.666667 0],'Marker','o', 'MarkerSize',6);
% xlabel('Data');
% ylabel('Probability')


% --- Create fit "LogNormFit"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('lognormal',[ -1.470473176121, 0.4677946640053])
pd1 = fitdist(arg_1,FitType); %'lognormal'
hLine = probplot(gca,pd1);
set(hLine,'LineStyle','-', 'LineWidth',2,'Color',col);%'Color',[1 0 0],

% --- Create fit "LognormFit"

% % Fit this distribution to get parameter values
% % To use parameter estimates from the original fit:
% %     pd2 = ProbDistUnivParam('lognormal',[ -1.351826985037, 0.4960531095293])
% pd2 = fitdist(arg_2, 'lognormal');
% hLine = probplot(gca,pd2);
% set(hLine,'Color',[0 0 1],'LineStyle','-', 'LineWidth',2);
PD = [PD,pd1];
end
% Adjust figure
box on;
grid on;
%hold off;
