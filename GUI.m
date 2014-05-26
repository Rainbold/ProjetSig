function varargout = GUI(varargin)
%GUI M-file for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('Property','Value',...) creates a new GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI('CALLBACK') and GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 18-May-2014 19:32:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Initialisation of some variables
handles.ecg_raw = [];
handles.ecg_data = [];
handles.ecg_Fs = 0;
handles.offset = 0;
handles.offset_freq = 0;
handles.time_window = 1440;
handles.freq_window = 4000;
handles.P = [];
handles.Q = [];
handles.R = [];
handles.S = [];
handles.T = [];
handles.method = 1;
handles.threshold_R = 40;
handles.threshold_QS = 1/12;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_plot_time(handles)
axes(handles.axes1);
N = length(handles.ecg_data); % N = seconds * Fs
cla(handles.axes1);
of = handles.offset;

handles.ecg_data = denoising(handles, handles.offset, handles.time_window);

if(get(handles.time_display_button,'Value') == 1)
    plot((of:(of+N-1))/handles.ecg_Fs, handles.ecg_data);
    hold on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Time evolution');
    if(N > 1)
        xlim([(of/handles.ecg_Fs) ((of+N-1)/handles.ecg_Fs)]);
    end
end

[ P, Q, R, S, T ] = update_pqrst(handles);
handles.P = P;
handles.Q = Q;
handles.R = R;
handles.S = S;
handles.T = T;

% P
if(get(handles.P_button,'Value') == 1)
    plot((P-1+of)/handles.ecg_Fs, handles.ecg_data(P),'+g');
    hold on;
end
% Q
if(get(handles.Q_button,'Value') == 1)
    plot((Q-1+of)/handles.ecg_Fs, handles.ecg_data(Q),'+c');
    hold on;
end
% R
if(get(handles.R_button,'Value') == 1)
    plot((R-1+of)/handles.ecg_Fs, handles.ecg_data(R),'+r');
    hold on;
end
% S
if(get(handles.S_button,'Value') == 1)
    plot((S-1+of)/handles.ecg_Fs, handles.ecg_data(S),'+m');
    hold on;
end
% T
if(get(handles.T_button,'Value') == 1)
    plot((T-1+of)/handles.ecg_Fs, handles.ecg_data(T),'+k');
    hold on;
end

crt = (1/(length(R)-1)) * sum(diff(R));
cardrhythm = 60 * handles.ecg_Fs / crt;

set(handles.txt_cardiac_rhythm, 'String', cardrhythm);

% Tachycardia
if(cardrhythm > 100)
   set(handles.text_tachy, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_tachy, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_tachy, 'String','Yes');
else
   set(handles.text_tachy, 'ForegroundColor',[0 0.5 0]); 
   set(handles.text_is_tachy, 'ForegroundColor',[0 0.5 0]);
   set(handles.text_is_tachy, 'String','No');
end

% Bradicardia
if(cardrhythm < 60)
   set(handles.text_brady, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_brady, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_brady, 'String','Yes');
else
   set(handles.text_brady, 'ForegroundColor',[0 0.5 0]); 
   set(handles.text_is_brady, 'ForegroundColor',[0 0.5 0]);
   set(handles.text_is_brady, 'String','No');
end

% Ectopic
dr = diff(R);
ddr = diff(dr);
threshold = 15;
if(ddr(ddr > threshold) ~= 0)
   set(handles.text_ectopic, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_ectopic, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_ectopic, 'String','Yes');
else
   set(handles.text_ectopic, 'ForegroundColor',[0 0.5 0]); 
   set(handles.text_is_ectopic, 'ForegroundColor',[0 0.5 0]);
   set(handles.text_is_ectopic, 'String','No');
end

% AF
if(length(P) < (1/10) * length(R))
   set(handles.text_AF, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_AF, 'ForegroundColor',[1 0 0]);
   set(handles.text_is_AF, 'String','Yes');
else
   set(handles.text_AF, 'ForegroundColor',[0 0.5 0]); 
   set(handles.text_is_AF, 'ForegroundColor',[0 0.5 0]);
   set(handles.text_is_AF, 'String','No');
end

update_respiratory_bpm(handles);

guidata(findobj('Tag', 'figure1'), handles);


function update_plot_freq(handles)
axes(handles.axes2);
of = handles.offset_freq;
N = handles.freq_window;
cla(handles.axes2);

data = denoising(handles, of, N);

if(get(handles.freq_display_button, 'Value') == 1);
    if(~isempty(handles.ecg_raw))
        dsp_ecg = abs(fft(data));
    else
        dsp_ecg = [];
    end
    plot(((0:(N-1))/N-0.5)*handles.ecg_Fs, fftshift(dsp_ecg));
    title('Power spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
end

function [ P, Q, R, S, T ] = update_pqrst(handles)
if(~isempty(handles.ecg_data))
    [ P, Q, R, S, T ] = PQRST(handles.method, handles.ecg_data, handles.ecg_Fs, handles.threshold_R, handles.threshold_QS);
else
    P = [];
    Q = [];
    R = [];
    S = [];
    T = [];
end

function ecg  = denoising(handles, offset, window)
if(~isempty(handles.ecg_raw))
    slp = get(handles.slider_lp, 'Value');
    shp = get(handles.slider_hp, 'Value');
    sRdb = get(handles.slider_Rdb, 'Value');         
    
    [bl, al] = butter(2, slp ^ 3, 'low');
    [bh, ah] = cheby1(2, shp ^ 3, sRdb ^ 3, 'high');
    
    data = handles.ecg_raw((1+offset):(offset+window));
    switch(get(handles.lp_button, 'Value') + get(handles.hp_button, 'Value') * 2)
        case 0
            ecg = data;
        case 1 % Low-pass filter
            ecg = filter(bl, al, data);
        case 2 % High-pass filter
            ecg = filter(bh, ah, data);
        case 3 % Low-pass and High-pass filter
            ecg = filter(bl, al, data);
            ecg = filter(bh, ah, ecg);
    end
end

function update_respiratory_bpm(handles)
R = handles.R;
if(length(R) >= 3)
    Fs = handles.ecg_Fs;

    dr = diff(R);

    N = length(dr);
    v = [];
    for i=1:N-1

        for j=R(i):R(i+1)
            v = [ v (dr(i) + ( (dr(i+1)-dr(i) ) .* (R(i+1) - R(i)).^-1 ) * (j-R(i))) ];
        end
    end
    N2 = 10^5;
    sp = fftshift(abs(fft((v-mean(v)),N2)));
    half = sp((N2/2+1):end);
    % N/2 -> Fs/2 Hz
    % X -> 0.4Hz
    % X = N/2 * 0.4 / (Fs/2)
    data = half(1:floor((N2/2) * 0.4 / (Fs/2)));
    d = diff(data);
    dd = diff(d);
    t1=d(1:end-1);
    t2=d(2:end);
    tt=t1.*t2;
    indx=find(tt<0);
    indx(dd(indx) > 0) = [];
    if(length(indx) >= 4)
        f = indx(4) * 60 * Fs / N2;
        set(handles.txt_respiratory_rhythm, 'String', f);
    else
        set(handles.txt_respiratory_rhythm, 'String', '- -');
    end
end

% --- Executes during object creation, after setting all properties.
function display_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider_offset_Callback(hObject, eventdata, handles)
% hObject    handle to slider_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.offset = floor(get(hObject, 'Value') * (length(handles.ecg_raw)-handles.time_window));
if(~isempty(handles.ecg_raw))
    handles.ecg_data = handles.ecg_raw((1+handles.offset):(handles.time_window+handles.offset));
end
guidata(hObject, handles);
update_plot_time(handles);



% --- Executes during object creation, after setting all properties.
function slider_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in P_button.
function P_button_Callback(hObject, eventdata, handles)
% hObject    handle to P_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_button

update_plot_time(handles);

% --- Executes on button press in Q_button.
function Q_button_Callback(hObject, eventdata, handles)
% hObject    handle to Q_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Q_button

update_plot_time(handles);

% --- Executes on button press in R_button.
function R_button_Callback(hObject, eventdata, handles)
% hObject    handle to R_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of R_button

update_plot_time(handles);


% --- Executes on button press in S_button.
function S_button_Callback(hObject, eventdata, handles)
% hObject    handle to S_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of S_button

update_plot_time(handles);

% --- Executes on button press in T_button.
function T_button_Callback(hObject, eventdata, handles)
% hObject    handle to T_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T_button

update_plot_time(handles);

% --- Executes on slider movement.
function slider_threshold_R_Callback(hObject, eventdata, handles)
% hObject    handle to slider_threshold_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.threshold_R = get(hObject, 'Value') * 100;
if(handles.threshold_R <= 0)
   handles.threshold_R = 1; 
end
guidata(hObject, handles);
update_plot_time(handles);



% --- Executes during object creation, after setting all properties.
function slider_threshold_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_threshold_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', 40/100);


% --- Executes on slider movement.
function slider_threshold_QS_Callback(hObject, eventdata, handles)
% hObject    handle to slider_threshold_QS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.threshold_QS = get(hObject, 'Value') * (1/2);

guidata(hObject, handles);
update_plot_time(handles);



% --- Executes during object creation, after setting all properties.
function slider_threshold_QS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_threshold_QS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', (1/12) * 2);


% --- Executes on button press in time_display_button.
function time_display_button_Callback(hObject, eventdata, handles)
% hObject    handle to time_display_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of time_display_button

update_plot_time(handles);


% --- Executes on button press in freq_display_button.
function freq_display_button_Callback(hObject, eventdata, handles)
% hObject    handle to freq_display_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of freq_display_button

update_plot_freq(handles);


% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.time_window = floor((get(hObject, 'Value')^2) * length(handles.ecg_raw));
if(handles.time_window == 0)
    handles.time_window = 2;
end
set(handles.edit_time_interval, 'String', handles.time_window/handles.ecg_Fs);

if(handles.offset+handles.time_window > length(handles.ecg_raw))
   handles.offset = length(handles.ecg_raw) - handles.time_window;
end
if(~isempty(handles.ecg_raw ~= 0))
    handles.ecg_data = handles.ecg_raw((1+handles.offset):(handles.offset+handles.time_window));
end

set(handles.slider_offset, 'Value', handles.offset / (length(handles.ecg_raw) - handles.time_window));

cla(handles.axes1, 'reset');
guidata(hObject, handles);
update_plot_time(handles);



% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_freq_Callback(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.freq_window = floor((get(hObject, 'Value')^2) * length(handles.ecg_raw));
if(handles.freq_window == 0)
    handles.freq_window = 2;
end

set(handles.edit_freq_interval, 'String', handles.freq_window/handles.ecg_Fs);

if(handles.offset_freq+handles.freq_window > length(handles.ecg_raw))
   handles.offset_freq = length(handles.ecg_raw) - handles.freq_window;
end

set(handles.slider_offset_freq, 'Value', handles.offset_freq / (length(handles.ecg_raw) - handles.freq_window));

cla(handles.axes2, 'reset');
guidata(hObject, handles);
update_plot_freq(handles);


% --- Executes during object creation, after setting all properties.
function slider_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function open_ecg_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open_ecg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile('*.mat');
if ~isequal(filename, 0)
    ecg_struct = load([filepath filename]);
   
    handles.ecg_raw = ecg_struct.ecg;
    handles.ecg_Fs = ecg_struct.Fs;
    handles.offset = 0;
    handles.time_window = 1440;
    handles.freq_window = 4000;
    handles.ecg_data = ecg_struct.ecg((1+handles.offset):(handles.offset+handles.time_window));

    guidata(hObject, handles);
    
    update_plot_time(handles);
    update_plot_freq(handles);
    
    set(handles.file_name,'String',filename)
    
    set(handles.slider_time', 'Value', sqrt(handles.time_window/length(handles.ecg_raw)));
    set(handles.edit_time_interval, 'String', handles.time_window/handles.ecg_Fs);
    
    set(handles.slider_freq', 'Value', sqrt(handles.freq_window/length(handles.ecg_raw)));
    set(handles.edit_freq_interval, 'String', handles.freq_window/handles.ecg_Fs);
    
    set(handles.slider_offset, 'Value', 0);
    set(handles.slider_offset_freq, 'Value', 0);
    
    set(handles.text_Fs, 'String', sprintf('%s%d', 'Fs : ', handles.ecg_Fs));
end



function edit_time_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_time_interval as a double

handles.time_window = floor(str2double(get(hObject,'String')) * handles.ecg_Fs);

set(handles.slider_time, 'Value', sqrt(handles.time_window/length(handles.ecg_raw)));

handles.ecg_data = handles.ecg_raw((1+handles.offset):handles.time_window);

guidata(hObject, handles);
update_plot_time(handles);


% --- Executes during object creation, after setting all properties.
function edit_time_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freq_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freq_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freq_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_freq_interval as a double

handles.freq = floor(str2double(get(hObject,'String')) * handles.ecg_Fs);

set(handles.slider_freq, 'Value', sqrt(handles.freq_window/length(handles.ecg_raw)));

guidata(hObject, handles);
update_plot_freq(handles);

% --- Executes during object creation, after setting all properties.
function edit_freq_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freq_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function time_display_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_display_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Value', 1);


% --- Executes during object creation, after setting all properties.
function method1_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Value', 1);


% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.NewValue
    case handles.method1_button
        handles.method = 1;
    case handles.method2_button
        handles.method = 2;
    case handles.method3_button
        handles.method = 3;
end
guidata(hObject, handles);
update_plot_time(handles);



% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes1, 'reset');
update_plot_time(handles);


% --- Executes on button press in clear_freq.
function clear_freq_Callback(hObject, eventdata, handles)
% hObject    handle to clear_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes2, 'reset');
update_plot_freq(handles);


% --- Executes on button press in lp_button.
function lp_button_Callback(hObject, eventdata, handles)
% hObject    handle to lp_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lp_button

update_plot_time(handles);
update_plot_freq(handles);


% --- Executes on button press in hp_button.
function hp_button_Callback(hObject, eventdata, handles)
% hObject    handle to hp_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hp_button

update_plot_time(handles);
update_plot_freq(handles);

% --- Executes on slider movement.
function slider_lp_Callback(hObject, eventdata, handles)
% hObject    handle to slider_lp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if(get(handles.lp_button, 'Value') == 1)
    update_plot_time(handles);
    update_plot_freq(handles);
end
f = (get(handles.slider_lp, 'Value')^2) * handles.ecg_Fs / 2;
set(handles.text_cp_lp, 'String', sprintf('%5.2f%s', f, ' Hz')); 


% --- Executes during object creation, after setting all properties.
function slider_lp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_lp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_hp_Callback(hObject, eventdata, handles)
% hObject    handle to slider_hp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if(get(handles.hp_button, 'Value') == 1)
    update_plot_time(handles);
    update_plot_freq(handles);
end
f = (get(handles.slider_hp, 'Value')^2) * handles.ecg_Fs / 2;
set(handles.text_cp_hp, 'String', sprintf('%5.2f%s', f, ' Hz')); 


% --- Executes during object creation, after setting all properties.
function slider_hp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_hp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_Rdb_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Rdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if(get(handles.hp_button, 'Value') == 1)
    update_plot_time(handles);
    update_plot_freq(handles);
end


% --- Executes during object creation, after setting all properties.
function slider_Rdb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Rdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_offset_freq_Callback(hObject, eventdata, handles)
% hObject    handle to slider_offset_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.offset_freq = floor(get(handles.slider_offset_freq, 'Value') * (length(handles.ecg_raw) - handles.freq_window));

guidata(hObject, handles);
update_plot_freq(handles);


% --- Executes during object creation, after setting all properties.
function slider_offset_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_offset_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in sync_button.
function sync_button_Callback(hObject, eventdata, handles)
% hObject    handle to sync_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider_freq, 'Value', get(handles.slider_time, 'Value'));
set(handles.edit_freq_interval, 'String', get(handles.edit_time_interval, 'String'));
set(handles.slider_offset_freq, 'Value', get(handles.slider_offset, 'Value'));

handles.offset_freq = handles.offset;
handles.freq_window = handles.time_window;

cla(handles.axes2, 'reset');
update_plot_freq(handles);
