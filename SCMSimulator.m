function varargout = SCMSimulator(varargin)
% SCMSIMULATOR MATLAB code for SCMSimulator.fig
%      SCMSIMULATOR, by itself, creates a new SCMSIMULATOR or raises the existing
%      singleton*.
%
%      H = SCMSIMULATOR returns the handle to a new SCMSIMULATOR or the handle to
%      the existing singleton*.
%
%      SCMSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCMSIMULATOR.M with the given input arguments.
%
%      SCMSIMULATOR('Property','Value',...) creates a new SCMSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCMSimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCMSimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCMSimulator

% Last Modified by GUIDE v2.5 30-Sep-2018 15:15:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCMSimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @SCMSimulator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% 不仅仅这部分改写会造成错误，而且 一般性的createFcn也会出错，重新生成了对象，无法找到原始对象；

% --- Executes just before SCMSimulator is made visible.
function SCMSimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCMSimulator (see VARARGIN)

% Choose default command line output for SCMSimulator
handles.output = hObject;

axes(handles.ax1);
I = imread('uts.jpg');
imshow(I);

axes(handles.ax2)
Inj = imread('njupt.jpg');
imshow(Inj);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using SCMSimulator.
%  if strcmp(get(hObject,'Visible'),'off')
%      plot(rand(5));
%  end



% UIWAIT makes SCMSimulator wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCMSimulator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
pop_Ind = get(handles.popupmenu1,'Value');
 set(0,'DefaultFigureVisible', 'on') 
 axes(handles.ax3);
 cla
 x = [1:300];
  CIR_Gen = handles.CIR_Gen;
switch(pop_Ind)
    case 1
   % axes(handles.ax3)
    plot(x*0.4,CIR_Gen(x,50));
    grid on;

    case 2
           plot(x*0.4,handles.CIR_Gen2(x,50));
    grid on; 
    case 3
            plot(x*0.4,handles.CIR_Gen3(x,50));
    grid on;
    case 4
      grid on;      plot(x*0.4,handles.CIR_Gen(x,50),'r'); hold on; grid on;
       plot(x*0.4,handles.CIR_Gen2(x,50),'b'); plot(x*0.4,handles.CIR_Gen3(x,50),'g'); hold off;  
       legend('SCM','SCM-PSC','SCM-FSC')
end
  xlim([0 120]);
  xlabel('Generated Channel Impulse Response (ns)');
  set(0,'DefaultFigureVisible', 'off') 
%   set(axes1,'FontSize',12,'YTick',...
%     [-3.71901648545568 -3.09023230616781 -2.32634787404084 -1.2815515655446 -0.674489750196082 0 0.674489750196082 1.2815515655446 2.32634787404084 3.09023230616781 3.71901648545571],...
%     'YTickLabel',...
%     {'0.0001','0.001','0.01','0.1','0.25','0.5','0.75','0.9','0.99','0.999','0.9999'});
  
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Plot SCM','Plot SCM-PSC','Plot SCM-FSC','Plot All'});
%{'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on mouse press over axes background.
function ax1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = imread('UTSLogo.png');
imshow(I);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% axes(handles.ax1);
% cla;
% popup_sel_index = get(handles.popupmenu1, 'Value');
% switch popup_sel_index
%     case 1
%         plot(rand(5));
%     case 2
%         plot(sin(1:0.01:25.99));
%     case 3
%         bar(1:.5:10);
%     case 4
%         plot(membrane);
%     case 5
%         surf(peaks);
% end
% I = imread('UTSLogo.png');
% imshow(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 尽可能简单进行SparseChannelModeling RMS分析 SMV_SCMSCA
% SMV-SPC 仿真程序，
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CopyRight@ %%%%%%%%%%%
% PengfeiCui
% 2018-09-29 UTS Building 11， 
% FEIT，School of communication & computation
% PerfeyCui@126.com
% 对应的论文：“Sparse Channel Modeling for Measured and Simulated Wireless Propagation
% Scenarios”，Peng-Fei Cui, J. Andrew Zhang, Wen-Jun Lu, Y. Jay Guo, Hong-Bo Zhu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Step 1 加载数据，计算RMS，利用TDL（直接采样20点，TDL平均20点）
%Step 2 统计某数据集的，稀疏特性，比例值分布，shape分布，非0元分布
        % 生成一个参数列表； 注意一个约束 sum=1；
%Step 3 仿真生成的信道模型，利用变稀疏矢量方式，三个模型；
        % 1）生成统计量，2）场景指数vec 统计量，3）稀疏计算量
        
%Step 4 对比实验，用TDL的随机生成的信号20taps和SCM对比；
        % 这个效果必然好太多吧！

  %% Final : a loop for recording RMS data;
  name={'OMP_wav_Dis','BP_wav_Dis','OMP_wav_Bonbo','BP_wav_bonbo','OMP_wav_Hei','BP_wav_Hei',...
      'OMP_wav_SV','BP_wav_SV','OMP_wav_Ray','BP_wav_Ray',};
        
  for iiAll =1%1:10 %3:10
      data_ind = get(handles.pm_data,'Value');
      alg_ind = get(handles.pm_algorithm,'Value');
      iiAll=2*(data_ind-1)+alg_ind
% Step1
switch(iiAll)
    case 1
%1 OMP Wavelet BanDis 400点数据
load('SCMAnalysis BanDis10x64_Wavlet_csOMP_201807Final_001.mat')
    case 2
%2 BP wave 
load('SCMAnalysis BanDis10x64_Wavlet_csBPm_201807Final_001.mat')
    case 3
%1b OMP wavelet BANBONBO;
load('SCMAnalysis BanBONBO10x80_Wavlet_csOMP_201807Final_001.mat')
    case 4
% 2b BP wave BAN bonbo
load('SCMAnalysis BanBONBO10x80_Wavlet_csBPm_201807Final_001_2.mat')
    case 5
%1c OMP wave BAN height
load('SCMAnalysis BanHeightVariation12x40_Wavlet_csOMP_201807Final_001.mat')
    case 6
%2c BP wave BAN Hei
load('SCMAnalysis BanHeightVariation12x40_Wavlet_csBPm_201807Final_001.mat')
    case 7
%3 sv modeling
load('SCMAnalysis SVgeneratedData_Wavlet_csOMP_201807Final_001.mat')
    case 8
load('SCMAnalysis SVgeneratedData_Wavlet_csBPm_201807Final_001.mat')
    case 9
%4 raymodeling: omp bp, waveletdic
load('SCMAnalysis ExpDecaying_Wavlet_csOMP_201807Final_001.mat')
    case 10
load('SCMAnalysis ExpDecaying_Wavlet_csBPm_201807Final_001.mat')

end

load('Sig_TDL.mat') % code 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 不显示图片；
  set(0,'DefaultFigureVisible', 'off') 
% Step 2 Statistics for sparse vector;
 dictionary1 = {{'sym3',5}};%{'RnIdent'}% %{dictionary{ii}};  % % ,'dct','RnIdent','poly',  ,'dct','poly','RnIdent' %{'sin'} ,{'poly'},{'RnIdent'} %{'sym4',5}
[mpdict,nbvect,LST,LONGS] = wmpdictionary(500,'lstcpt',dictionary1); 
LevIndR = LONGS{1, 1}(1:end-1,1);
LevInd = cumsum(LevIndR)

%I-3 分析小波树3（全部tap分析) ： 分析出tap1的 1）levl号码； 2）数值；
% FileName = {'C:\Users\perfe\Documents\UTS_Work\2017JulAugAu\20170910\ImpulseWaveletMatlab\SparseChannAnalysisResults\BanSparseChanAnalysis BanBONBO10x80SCANmpcFilter_201807Test1sym35.mat'}
% %{'C:\Users\perfe\Documents\UTS_Work\2017JulAugAu\20170910\ImpulseWaveletMatlab\SparseChannAnalysisResults\BanSparseChanAnalysis BanBONBO10x80SCANmpcFilter_201806Test.mat'};
% %[pathstr,name,ext] = fileparts(filename) 
%  % 耗计算，单独列， tackleSparseCoeff_Pro.mlx  XXX
% load(FileName{1})
[m,n]=size(SparsityChAnalysisResults);

Nmin = min(cell2mat(SparsityChAnalysisResults(:,5)));
Nmax = max(cell2mat(SparsityChAnalysisResults(:,5)));

Tap1Val = zeros(500,m); Tap1LevRaw = zeros(500,m);
Tap1Lev=zeros(500,m);  Tap1Loc = zeros(500,m);
for jjj = 1:400%:m %400
   %[ValT,IndT] = max(abs(SparsityChAnalysisResults{jjj, 22}) ); 
   [ValT,IndT] = sort(abs(SparsityChAnalysisResults{jjj, 22}),'descend' ); 
   InT = find(ValT==0);%InT(1)为实际的tap数目；
   
   Tap1Val(1:InT(1),jjj) = ValT(1:InT(1));%[Tap1Val;ValT];
   Tap1LevRaw((1:InT(1)),jjj) = IndT(1:InT(1));%[];
   % 只有Sym4-5是隔断点矢量：[17 34 66 128 253 502]
  [Tap1Loc(1:InT(1),jjj), Tap1Lev(1:InT(1),jjj) ] = SubFun_judgeLev...
      (Tap1LevRaw(1:InT(1),jjj),LevInd');
end
% %Draw
% for ii = 1:20;%Nmax
%  figure,  hold on;
%  histogram(Tap1Lev(ii,:))
% legend(num2str(ii)) 
% % end
%  'success'

% Step2 b, 解析shape值（level的）构成： 倒置的1-6 变为 1 6 5 4 3 2
%生成时候需要反过来，倒转2-6生成；
% Inverse the sequence;
TapLevInv = 8-Tap1Lev; %（用7减）倒置，0-》7,6-》1；
xx = find(TapLevInv==7); %误统计的7，其实是0，无所谓的；先去掉吧；
TapLevInv(xx) =1;%0.1;

PD_LevInv_save = []; PD_LevInv_saveSigma = [];
for jj = 1:20
    figure,t = TapLevInv(jj,:);
    xx = find(t<=7);
    temp = t(xx);
   [pd1,pd2,pd3] = createFitOfTap1to3(temp,'Noo'); 
   PD_LevInv_save=[PD_LevInv_save,pd1];
   PD_LevInv_saveSigma=[PD_LevInv_saveSigma,pd1.sigma];
   xlabel(['Data ' num2str(jj)])
end


%%%%% location lognorm analysis%%%%%%%%%%%%%%%%%%%%%%%%
TapLoc = Tap1Loc;
xx = find(Tap1Loc==0);
%TapLoc(xx)=1e-4;

PD_Loc_save = []; PD_Loc_saveSigma = [];PD_Loc_savemu = [];
for jj = 1:20
    figure %TapLevInv,
   %[pd1,pd2,pd3] = createFitOfTap1to3(TapLevInv(jj,:)) 
%    [pd1] = createFit1to1(TapLevInv(jj,:),'norm','cdf');
  t = TapLoc(jj,:);
    xx = find(t<1);
    temp = t(xx); xx2 = find(temp>0);
    temp2 = temp(xx2);
[pd1] = createFit1to1(temp2,'lognorm','pro'); %log
   %TapLocPDsave=[TapLocPDsave,pd1];
    PD_Loc_save=[PD_Loc_save,pd1];
   PD_Loc_savemu=[PD_Loc_savemu,pd1.mu];
   PD_Loc_saveSigma=[PD_Loc_saveSigma,pd1.sigma];
   xlabel(['Data ' num2str(jj)])
end

%%%%% tap values：
PD_TapVal_save = []; %除了前1-5个，基本符合lognorm分布；
PD_Tap_saveSigma = [];PD_Tap_savemu = [];
TapVal = Tap1Val;
xx = find(Tap1Val==0);
TapVal(xx)=1e-4;
for jj = 1:20
    TapTem = TapVal(jj,:);
    xx = find(TapTem>1e-3);
    TapTem = TapTem(xx);
    figure %TapLevInv,
    [pd1] = createFit1to1(TapTem,'norm','pro'); %log
   PD_TapVal_save=[PD_TapVal_save,pd1];
   PD_Tap_saveSigma = [PD_Tap_saveSigma,pd1.sigma];
   PD_Tap_savemu = [PD_Tap_savemu,pd1.mu];
   xlabel(['Data ' num2str(jj)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% table:
  set(0,'DefaultFigureVisible', 'on') 
  t = handles.uitable;
  %m = magic(3);
  set(t,'Data',[PD_Tap_savemu',PD_Tap_saveSigma',PD_Loc_savemu',PD_Loc_saveSigma',PD_LevInv_saveSigma']);
  set(t,'ColumnName',{'Coeff\mu';'Coef\sig';'Loc\mu';'Loc\Sig';'Shape\Sig'})
  %{'CoefNorm\mu';'Coef\sig';'ShapeHnorm \Sigma';'LocLogNorm\mu';'LocLogNorm\Sigma'})
  
    set(0,'DefaultFigureVisible', 'off') 
%Step 3 信道生成模型： %% Generation CIR
%% Generation CIR
CoeffGen = [];CoeffIndGen = []; LevelSeqGen=[]; GenNum = 400;
 for pdi = 1:20
Coeff_Gen(pdi,:)=random(PD_TapVal_save(pdi),1,GenNum);
Coeff_Loc_Gen(pdi,:)=random(PD_Loc_save(pdi),1,GenNum); 
Coeff_LevInv_Gen(pdi,:)=random(PD_LevInv_save(pdi),1,GenNum); 
%每一行是一种参数值，如1st是第一个tab的值；每一列是一个系数组；
 end

 % tackle coeff numbers：
 xxx = find(Coeff_Loc_Gen>1); size(xxx);
 Coeff_Loc_Gen(xxx) = 1; clear xxx;
 
 % 1)取整、限界 2）顺序倒置，除了1以外， 2-6,6-》2
% Coeff_LevInv_Gen=floor(Coeff_LevInv_Gen);
Coeff_LevInv_Gen=round(Coeff_LevInv_Gen);
xxx = find(Coeff_LevInv_Gen<1);size(xxx);
Coeff_LevInv_Gen(xxx) = 1; clear xxx;
xxx = find(Coeff_LevInv_Gen>6);size(xxx);
Coeff_LevInv_Gen(xxx) = 6; clear xxx;

xxx = find(Coeff_LevInv_Gen>1);size(xxx)
Coeff_LevInv_Gen(xxx) = 8-Coeff_LevInv_Gen(xxx); clear xxx;
% 统计规律调整，tap1 tap2全部为scale1； %%%%%
Coeff_LevInv_Gen(1:2,:) = 1;
%Coeff_LevInv_Gen(1,:) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Method 2: coeffGen2 ae^-b+c;
x=[1:20];
y = 0.78*exp(-0.35*x)+0.02;
Coeff_Gen2=repmat(y',1,400);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Method 3: coeffGen2 ae^-b+c;
x3=[1:20];

y1 = Coeff_Gen(1,1);
yk = Coeff_Gen(20,1);
cc = yk;
bb = log((yk*(y1-yk)+sqrt(yk^2*(y1-yk)^2+(1-20*yk)^2*(1-y1^2-(20-1)*yk^2)))/(1-20*yk^2)^2);
aa = (y1-yk)/(exp(bb));
y3 = aa*exp(bb*x3)+cc;

Coeff_Gen3=repmat(y3',1,400);


%Step 3:系数生成：
Coeff_Ind = zeros(20,GenNum); 
% 生成系数值，生成起始值，
xxx = find(Coeff_LevInv_Gen>1);
Coeff_Ind(xxx)=LevInd(Coeff_LevInv_Gen(xxx)-1); clear xxx;
Coeff_Ind=Coeff_Ind+LevIndR(Coeff_LevInv_Gen).* Coeff_Loc_Gen;
% 取整后，排除异常值；
Coeff_Ind2 = round(Coeff_Ind);
xxx = find(Coeff_Ind2<1);
Coeff_Ind2(xxx) = 1; clear xxx;
xxx = find(Coeff_Ind2>502);
Coeff_Ind2(xxx) = 502; clear xxx;

% Gen
CIR_Gen = zeros(500,GenNum);
for jjj = 1:GenNum
CIR_Gen(:,jjj) =(mpdict(:,Coeff_Ind2(:,jjj))* Coeff_Gen(:,jjj) );
end
figure,plot(CIR_Gen(:,50))
figure,mesh(CIR_Gen)

view([-79.594 5.314])
figure,mesh(Coeff_Gen)

% Gen2
CIR_Gen2 = zeros(500,GenNum);
for jjj = 1:GenNum
CIR_Gen2(:,jjj) =(mpdict(:,Coeff_Ind2(:,jjj))* Coeff_Gen2(:,jjj) );
end
figure,plot(CIR_Gen2(:,50))
figure,mesh(CIR_Gen2)
figure,mesh(Coeff_Gen2)

% Gen3
CIR_Gen3 = zeros(500,GenNum);
for jjj = 1:GenNum
CIR_Gen3(:,jjj) =(mpdict(:,Coeff_Ind2(:,jjj))* Coeff_Gen3(:,jjj) );
end
figure,plot(CIR_Gen3(:,50))
figure,mesh(CIR_Gen3)
figure,mesh(Coeff_Gen3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.CIR_Gen = CIR_Gen;
handles.CIR_Gen2 = CIR_Gen2;
handles.CIR_Gen3 = CIR_Gen3;
guidata(hObject,handles)
  set(0,'DefaultFigureVisible', 'on') 
axes(handles.ax3)
x = [1:300];
plot(x*0.4,CIR_Gen3(x,50));
grid on; 
xlabel('Generated Channel Impulse Response (ns)');
%xlabel('Generated CIRs');
xlim([0 120])
%ylim([-0.05 0.4])
  set(0,'DefaultFigureVisible', 'off') 

%Step 4: check; KPI 参量怎么样?
RMS_Rec = []; % 稀疏恢复信号的RMS值，分布
D_Rec = [];
for ii = 1:400
rec_sig = SparsityChAnalysisResults{ii, 20};
timeResolution = 0.4; %ns
[rms_rec,dl] = RMS_DS(rec_sig,timeResolution);
RMS_Rec = [RMS_Rec,rms_rec];
D_Rec = [D_Rec,dl];
end

RMS_tdl = []; % TDL等间隔的信号的RMS值，分布
D_tdl =[];  
Sig_tdl= Sig_TDL(1+(iiAll-1)*20:iiAll*20,:);%[];
for iij = 1:400
rec_sig_tdl=Sig_tdl(:,iij);
timeResolution = 0.4*25; %ns
[rms_rec,dl] = RMS_DS(rec_sig_tdl,timeResolution);
RMS_tdl = [RMS_tdl,rms_rec];
D_tdl = [D_tdl,dl];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCM modeling data
RMS_scm = []; % SCM等间隔的信号的RMS值，分布
D_scm =[];
for iij = 1:400
rec_sig = CIR_Gen(1:400,iij);%SparsityChAnalysisResults{iij, 20};
timeResolution = 0.4;%*25; %ns
[rms_rec,dl] = RMS_DS(rec_sig,timeResolution);
RMS_scm = [RMS_scm,rms_rec];
D_scm = [D_scm,dl];
end
figure,plot(RMS_Rec,RMS_scm,'r*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCM2 modeling data
RMS_scm2 = []; % SCM等间隔的信号的RMS值，分布
D_scm2 =[];
for iij = 1:400
rec_sig = CIR_Gen2(1:400,iij);%SparsityChAnalysisResults{iij, 20};
timeResolution = 0.4;%*25; %ns
[rms_rec,dl] = RMS_DS(rec_sig,timeResolution);
RMS_scm2 = [RMS_scm2,rms_rec];
D_scm2 = [D_scm2,dl];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCM3 modeling data
RMS_scm3 = []; % SCM等间隔的信号的RMS值，分布
D_scm3 =[];
for iij = 1:400
rec_sig = CIR_Gen3(1:450,iij);%SparsityChAnalysisResults{iij, 20};
timeResolution = 0.4;%*25; %ns
[rms_rec,dl] = RMS_DS(rec_sig,timeResolution);
RMS_scm3 = [RMS_scm3,rms_rec];
D_scm3 = [D_scm3,dl];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 画图 %%%%%%%%%%%%%%%%%%%%
figure,plot(RMS_Rec,RMS_tdl,'r*');
% 画RMS对比的 图；
figure,hold on, Ltxt=[];
probplot('normal'); % create empty plot of desired type

%for ii = 0%:2 % 0 means Rayleigh, 1 Rice 5, 2 Nakagami w=5,
    [PD] = createFitProbability1Fig(RMS_Rec,'lognorm','Yes','r');
   % PD.mu,PD.sigma, 
    Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
    
 [PD] = createFitProbability1Fig(RMS_scm','lognorm','Yes','b');
    %PD.mu,PD.sigma, 
 Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
  
  [PD] = createFitProbability1Fig(RMS_scm2','lognorm','Yes','k');
    %PD.mu,PD.sigma, 
 Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
 
  [PD] = createFitProbability1Fig(RMS_scm3','lognorm','Yes','y');
    %PD.mu,PD.sigma, 
 Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
 
 [PD] = createFitProbability1Fig(RMS_tdl','lognorm','Yes','g');
    %PD.mu,PD.sigma, 
 Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
 
 %%%%%%%%%%%%%%%%%%% temporialy save %%%%%%%%%%%%%
  % nam = [nam,'L1_Wavelet_banDis_FixVal'],rms = [rms;RMS_scm2];
  %%%nam = ['L1_Wavelet_banDis'],rms = [RMS_scm];
  % save('RMS_Gen201809in.mat','nam','rms')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end                  %Predicted SparseCoefficient PSC, Fitted SparseCoeff
xlabel(Ltxt),legend('Raw RMS','Lognorm, \mu=2.6,\sigma=0.5','SCM','Lognorm, \mu=2.8,\sigma=0.3','SCM-PSC','Lognorm, \mu=2.7,\sigma=0.3',...
    'SCM-FSC','Lognorm, \mu=2.9,\sigma=0.3','TDL','Lognorm, \mu=3.4,\sigma=0.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  set(0,'DefaultFigureVisible', 'on') 
axes(handles.ax4)
% 画Delay对比的 图；%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, Ltxt=[];  rmsPD = [];
probplot('normal'); % create empty plot of desired type

%for ii = 0%:2 % 0 means Rayleigh, 1 Rice 5, 2 Nakagami w=5,
    [PD] = createFitProbability1Fig(D_Rec,'lognorm','Yes','r');
    rmsPD = [rmsPD,PD.mu,PD.sigma];%PD.mu,PD.sigma, 
    Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
        [PD] = createFitProbability1Fig(D_scm,'lognorm','Yes','b');
    rmsPD = [rmsPD,PD.mu,PD.sigma];%PD.mu,PD.sigma, 
    Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
    [PD] = createFitProbability1Fig(D_scm2,'lognorm','Yes','k');
    rmsPD = [rmsPD,PD.mu,PD.sigma];%PD.mu,PD.sigma, 
    Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
    [PD] = createFitProbability1Fig(D_scm3,'lognorm','Yes','y');
    rmsPD = [rmsPD,PD.mu,PD.sigma];%PD.mu,PD.sigma, 
    Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
    
   [PD] = createFitProbability1Fig(D_tdl','lognorm','Yes','g');
   rmsPD = [rmsPD,PD.mu,PD.sigma]; %PD.mu,PD.sigma, 
 Ltxt = [Ltxt ' ' num2str(PD.mu) ' ' num2str(PD.sigma)];
%end
%xlabel(Ltxt),
xlabel('Root Mean Square Delay Spread (ns)')
%legend('Original data','SCM','TDL');
legend('Raw RMS','Lognorm','SCM','Lognorm','SCM-PSC','Lognorm',...
    'SCM-FSC','Lognorm','STDL','Lognorm');
set(handles.ax4,'FontSize',9,'YTick',...
    [-3.71901648545568 -3.09023230616781 -2.32634787404084 -1.2815515655446 -0.674489750196082 0 0.674489750196082 1.2815515655446 2.32634787404084 3.09023230616781 3.71901648545571],...
    'YTickLabel',...
    {'0.0001','0.001','0.01','0.1','0.25','0.5','0.75','0.9','0.99','0.999','0.9999'});
title('');
xlim([10 100]);
  set(0,'DefaultFigureVisible', 'off') 




'success'

 %%%%%%%%%%%%%%%%%%% statistical RMS save %%%%%%%%%%%%%
  % nam = [nam,'L1_Wavelet_banDis_FixVal'],rms = [rms;RMS_scm2];
  %nam = [1,{'OMP_Wavelet_banDis'}],rms = [RMS_Rec',RMS_scm',RMS_scm2',RMS_scm3',RMS_tdl'];
%  RMSPD = rmsPD;
% 
% nam = [nam;iiAll,name(iiAll)]; rms = [rms,RMS_Rec',RMS_scm',RMS_scm2',RMS_scm3',RMS_tdl'];
% RMSPD = [RMSPD;rmsPD];
%  save('RMS_GenPD201809.mat','nam','rms','RMSPD')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  set(0,'DefaultFigureVisible', 'on') 

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_sparsity.
function edit_sparsity_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_sparsity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_sparsity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sparsity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_measurement_Callback(hObject, eventdata, handles)
% hObject    handle to edit_measurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_measurement as text
%        str2double(get(hObject,'String')) returns contents of edit_measurement as a double


% --- Executes during object creation, after setting all properties.
function edit_measurement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_measurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_data.
function pm_data_Callback(hObject, eventdata, handles)
% hObject    handle to pm_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_data


% --- Executes during object creation, after setting all properties.
function pm_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pm_algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_algorithm.
function pm_algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to pm_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_algorithm


% --- Executes during object creation, after setting all properties.
function pm_dictionary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_dictionary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_dictionary.
function pm_dictionary_Callback(hObject, eventdata, handles)
% hObject    handle to pm_dictionary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_dictionary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_dictionary


% --- Executes on selection change in pm_scheme.
function pm_scheme_Callback(hObject, eventdata, handles)
% hObject    handle to pm_scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_scheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_scheme


% --- Executes during object creation, after setting all properties.
function pm_scheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pm_parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_parameter.
function pm_parameter_Callback(hObject, eventdata, handles)
% hObject    handle to pm_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_parameter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_parameter


% --- Executes when entered data in editable cell(s) in uitable.
function uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function editSimulationTimes_Callback(hObject, eventdata, handles)
% hObject    handle to editSimulationTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSimulationTimes as text
%        str2double(get(hObject,'String')) returns contents of editSimulationTimes as a double


% --- Executes during object creation, after setting all properties.
function editSimulationTimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSimulationTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mse_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mse as text
%        str2double(get(hObject,'String')) returns contents of edit_mse as a double


% --- Executes during object creation, after setting all properties.
function edit_mse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_sparsity_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_sparsity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_sparsity as text
%        str2double(get(hObject,'String')) returns contents of edit_max_sparsity as a double


% --- Executes during object creation, after setting all properties.
function edit_max_sparsity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_sparsity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
