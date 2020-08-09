function varargout = alzdetection_v2(varargin)
% ALZDETECTION_V2 MATLAB code for alzdetection_v2.fig
%      ALZDETECTION_V2, by itself, creates a new ALZDETECTION_V2 or raises the existing
%      singleton*.
%
%      H = ALZDETECTION_V2 returns the handle to a new ALZDETECTION_V2 or the handle to
%      the existing singleton*.
%
%      ALZDETECTION_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALZDETECTION_V2.M with the given input arguments.
%
%      ALZDETECTION_V2('Property','Value',...) creates a new ALZDETECTION_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alzdetection_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alzdetection_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alzdetection_v2

% Last Modified by GUIDE v2.5 23-Jun-2020 06:20:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alzdetection_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @alzdetection_v2_OutputFcn, ...
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


% --- Executes just before alzdetection_v2 is made visible.
function alzdetection_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alzdetection_v2 (see VARARGIN)

% Choose default command line output for alzdetection_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alzdetection_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alzdetection_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.infolist,'str',''); %Clears the info list box.
set(handles.plottingpanel,'visible','off');
handles.infolist.Max = length(get(handles.infolist,'str'));
handles.infolist.Value = handles.infolist.Max;

[filename,pathname] = uigetfile([pwd '\Dataset\*.mat'], 'Load Single Subject .mat file');
load(fullfile(pathname,filename));
set(handles.pathbox,'str',[pathname filename]); %Updates the path of the file into the pathbox
oldstr = get(handles.infolist,'str'); % The string as it is now.
addstr1 = {'------- Loading Data -------'};
addstr2 = {[strrep(filename,'.mat','') ' data loaded...']}; % The string to add to the stack.
set(handles.infolist,'str',[oldstr(:)',addstr1(:)',addstr2(:)']);  % Put the new string on top
handles.infolist.Max=length(get(handles.infolist,'str'));
handles.infolist.Value = handles.infolist.Max;

%% Plotting the selected data
% time_128 = 0:0.0078125:60-0.0078125;
newYlabels = {"F3","F4","C3","C4","P3","P4","O1","O2","F7","F8","T3","T4","T5","T6","Cz","Pz"};
% stackedplot(time_128,data,'Title',['Subject ',num2str(filename),' (Raw Data)'],'DisplayLabels',newYlabels);
% set(handles.plottingpanel,'visible','on');
data=data'; %update 'data' variable from (7680x16->16x7680)
map = vega20(100);
maxI = 0;
t = (0:size(data,2)-1)/128;% Fs = 128Hz
for i = 1 : size(data,1)
    normData = normalize_range(data(i,:), 0+(i-1), i);
    plot(handles.plotaxes,t,normData, 'linewidth',0.5,'color',map(28+i,:))
    hold(handles.plotaxes,'on');
    maxI = maxI+max(data(i,:));
end
hold(handles.plotaxes,'off');
axis tight
set(handles.plottingpanel,'visible','on');

%% Update Info List
oldstr = get(handles.infolist,'str');
addstr1 = {'--------- Data Info --------'};
addstr2 = {['Channels  : ' num2str(size(data,1))]};
addstr3 = {['Duration    : ' num2str(size(data,2))]};
addstr4 = {['------------------------------']};
set(handles.infolist,'str',[oldstr(:)',addstr1(:)',addstr2(:)',addstr3(:)',addstr4(:)']);
handles.infolist.Max=length(get(handles.infolist,'str'));
handles.infolist.Value = handles.infolist.Max;

handles.eeg_data=data'; %save after changing 'data' variable from (16x7680->7680x16)
guidata(hObject, handles);



function pathbox_Callback(hObject, eventdata, handles)
% hObject    handle to pathbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathbox as text
%        str2double(get(hObject,'String')) returns contents of pathbox as a double


% --- Executes during object creation, after setting all properties.
function pathbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in varbutton.
function varbutton_Callback(hObject, eventdata, handles)
% hObject    handle to varbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'eeg_data')
    msgbox('Load EEG before extracting VAR connectivity...');
    return;
end
set(handles.plottingpanel,'visible','on')
data = handles.eeg_data;
p  = 5; %VAR order (fixed)
ch = size(data,2);
t  = size(data,1);

oldstr = get(handles.infolist,'string');
addstr = {'Running VAR of Order 5 ....'};
set(handles.infolist,'str',[oldstr(:)',addstr(:)'])
handles.infolist.Max=length(get(handles.infolist,'string'));
handles.infolist.Value = handles.infolist.Max;

[VARtemp,Rr_Rrtemp] = varfit(p,data');
VARdata = reshape(VARtemp,ch,ch,p);

map = vega20(100);
delete(get(handles.plottingpanel,'child'))

for i = 1 : p
    h(i)=subplot(2,3,i,'Parent',handles.plottingpanel);%axis square;
    cla(h(i));
end
TickCh=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
LabelCh=['F3'; 'F4'; 'C3'; 'C4'; 'P3'; 'P4'; 'O1'; 'O2'; 'F7'; 'F8'; 'T3'; 'T4'; 'T5'; 'T6'; 'Cz'; 'Pz'];
for i = 1 : p
    imagesc(h(i),squeeze(VARdata(:,:,i)));%axis square;
    colormap(jet)
    %colorbar(h(i),'north');
    %caxis([-1 2]);
    
    h(i).XTick = TickCh;
    h(i).XTickLabel = LabelCh;
    h(i).YTick = TickCh;
    h(i).YTickLabel = LabelCh;
    title(h(i),['VAR Lag - ' num2str(i)]);
    xlabel(h(i),'From');ylabel(h(i),'To');
    set(h(i),'FontSize',8)
end
oldstr = get(handles.infolist,'string');
addstr1 = {'VAR coeffs estimation is done ....'};
addstr2 = {'------------------------------'};
set(handles.infolist,'str',[oldstr(:)',addstr1(:)',addstr2(:)'])
handles.infolist.Max=length(get(handles.infolist,'string'));
handles.infolist.Value = handles.infolist.Max;

% [d1,d2,d3] = size(VARdata); %VARdata dimension: 16x16x5
% VARdata = reshape(VARdata,d1,d2*d3); %reshape 16x16x5 -> 16x80

var_data.vartemp=VARtemp;
var_data.Rr = Rr_Rrtemp;
var_data.data=VARdata;
handles.var_data = var_data;
guidata(hObject, handles);



% --- Executes on button press in pdcbutton.
function pdcbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pdcbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'eeg_data')
    msgbox('Load EEG then click VAR before extracting PDC connectivity...');
    return;
end
if ~isfield(handles,'var_data')
    msgbox('Click VAR before extracting PDC connectivity...');
    return;
end
VARtemp = handles.var_data.vartemp;
Rr = handles.var_data.Rr;
p  = 5; %VAR order (fixed)
fs = 128; % sampling frequency
nfft=fs/2;
fc=fs/2;
ch=16;
fr = 1;
PDCdata = zeros(ch,ch,p);
[Ave_PDC,delta_PDC,theta_PDC,alpha_PDC,beta_PDC,gamma_PDC,gpdc,pdc,dtf,dc,coh,pcoh,pcoh2]=FunctionPDC2(VARtemp,Rr,nfft,fs,ch,p,fr);
% Store PDCdata
PDCdata(:,:,1) = delta_PDC; 
PDCdata(:,:,2) = theta_PDC;
PDCdata(:,:,3) = alpha_PDC;
PDCdata(:,:,4) = beta_PDC;
PDCdata(:,:,5) = gamma_PDC;

map = vega20(100);
delete(get(handles.plottingpanel,'child'))

for i = 1 : p
    h(i)=subplot(2,3,i,'Parent',handles.plottingpanel);%axis square;
    cla(h(i));
end
TickCh=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
LabelCh=['F3'; 'F4'; 'C3'; 'C4'; 'P3'; 'P4'; 'O1'; 'O2'; 'F7'; 'F8'; 'T3'; 'T4'; 'T5'; 'T6'; 'Cz'; 'Pz'];
% bands=['Delta';'Theta';'Alpha';'Beta ';'Gamma'];
bands=["Delta","Theta","Alpha","Beta","Gamma"];
for i = 1 : p
    imagesc(h(i),squeeze(PDCdata(:,:,i)));%axis square;
    colormap(jet)
    %colorbar(h(i),'north');
    %caxis([-1 2]);
    
    h(i).XTick = TickCh;
    h(i).XTickLabel = LabelCh;
    h(i).YTick = TickCh;
    h(i).YTickLabel = LabelCh;
%     title(h(i),['PDC band - ' bands(i)]);
    title(h(i),["PDC Band " + bands{i}]);
    xlabel(h(i),'From');ylabel(h(i),'To');
    set(h(i),'FontSize',8)
end

oldstr = get(handles.infolist,'string');
addstr1 = {'PDC coeffs estimation is done ....'};
addstr2 = {'------------------------------'};
set(handles.infolist,'str',[oldstr(:)',addstr1(:)',addstr2(:)'])
handles.infolist.Max=length(get(handles.infolist,'string'));
handles.infolist.Value = handles.infolist.Max;

pdc_data.data = PDCdata;
handles.pdc_data = pdc_data;
guidata(hObject, handles);



% --- Executes on selection change in infolist.
function infolist_Callback(hObject, eventdata, handles)
% hObject    handle to infolist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns infolist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from infolist


% --- Executes during object creation, after setting all properties.
function infolist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to infolist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadcnnbutton.
function loadcnnbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadcnnbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'var_data')
    msgbox('Load EEG and extract VAR first, at least...');
    return;
end
[filename,pathname] = uigetfile([pwd '\Models and Networks\*.hdf5'],'Choose CNN model for VAR/PDC?');
 model = importKerasNetwork(fullfile(pathname,filename));
 %model=load(fullfile(pathname,filename));
 %net=importKerasNetwork(model);
 
if contains(filename,'PDC') % PDC_CNN...
    if ~isfield(handles,'pdc_data')
        msgbox('Wrong model selection! Cannot classify PDC without Extracting PDC Connectivity...');
        return;
    end
end
 
 oldstr = get(handles.infolist,'string');
 addstr = {'Loading CNN Model ....'};
 addstr1 = {'DONE ....'};
 addstr2 = {'----------------------------'};
 set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr1(:)',addstr2(:)'])
 handles.infolist.Max=length(get(handles.infolist,'string'));
 handles.infolist.Value = handles.infolist.Max;
 
 model_data.model=model;
 model_data.filename = filename;
 handles.model_data = model_data;
guidata(hObject, handles);



% --- Executes on button press in cnnclassifybutton.
function cnnclassifybutton_Callback(hObject, eventdata, handles)
% hObject    handle to cnnclassifybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'model_data')
    msgbox('Load CNN Model first before classification...');
    return;
end
model = handles.model_data.model;
filename=handles.model_data.filename;
if contains(filename,'PDC') % PDC_CNN...
    pdc=handles.pdc_data.data;
    YPred = classify(model,pdc);
    
    if categorical (YPred) == '2'
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Non-Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    else
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    end
    
elseif contains(filename,'VAR') %VAR_CNN
    var=handles.var_data.data;
    YPred = classify(model,var);
    
    if categorical (YPred) == '2'
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Non-Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    else
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    end
end



% --- Executes on button press in loadannbutton.
function loadannbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadannbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'var_data')
    msgbox('Load EEG and extract VAR first, at least...');
    return;
end
[filename,pathname] = uigetfile([pwd '\Models and Networks\*.mat'],'Choose ANN network for VAR/PDC?');
load(fullfile(pathname,filename));

if contains(filename,'PDC') % PDC_CNN...
    if ~isfield(handles,'pdc_data')
        msgbox('Wrong model selection! Cannot classify PDC without Extracting PDC Connectivity...');
        return;
    end
end

 oldstr = get(handles.infolist,'string');
 addstr = {'Loading ANN Network ....'};
 addstr1 = {'DONE ....'};
 addstr2 = {'----------------------------'};
 set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr1(:)',addstr2(:)'])
 handles.infolist.Max=length(get(handles.infolist,'string'));
 handles.infolist.Value = handles.infolist.Max;
 
 net_data.net=net;
 net_data.filename = filename;
 handles.net_data = net_data;
 guidata(hObject, handles);



% --- Executes on button press in annclassifybutton.
function annclassifybutton_Callback(hObject, eventdata, handles)
% hObject    handle to annclassifybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'net_data')
    msgbox('Load ANN Network first before classification...');
    return;
end
net=handles.net_data.net;
filename=handles.net_data.filename;
if contains(filename,'PDC') % PDC_ANN...
    pdc=handles.pdc_data.data;
    pdc=pdc(:);
    YPred = net(pdc);
    
    if round(YPred) == 1
        oldstr = get(handles.infolist,'string');
        addstr = {'Status is : Non-Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    else
        oldstr = get(handles.infolist,'string');
        addstr = {'Status is : Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    end
    
elseif contains(filename,'VAR') %VAR_ANN
    var=handles.var_data.data;
    var=var(:);
    YPred = net(var);
    
    if round(YPred) == 1
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Non-Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    else
        oldstr = get(handles.infolist,'string');
        addstr = {'EEG Classified as : Alzheimer''s'};
        addstr2 = {'----------------------------'};
        addstr3 = {'RESTART APP FOR NEXT SUBJ.'};
        addstr4 = {'----------------------------'};
        set(handles.infolist,'str',[oldstr(:)',addstr(:)',addstr2(:)',addstr3(:)',addstr4(:)'])
        handles.infolist.Max=length(get(handles.infolist,'string'));
        handles.infolist.Value = handles.infolist.Max;
    end
end
