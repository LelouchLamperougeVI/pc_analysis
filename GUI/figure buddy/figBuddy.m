function varargout = figBuddy(varargin)
% FIGBUDDY MATLAB code for figBuddy.fig
%      FIGBUDDY, by itself, creates a new FIGBUDDY or raises the existing
%      singleton*.
%
%      H = FIGBUDDY returns the handle to a new FIGBUDDY or the handle to
%      the existing singleton*.
%
%      FIGBUDDY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGBUDDY.M with the given input arguments.
%
%      FIGBUDDY('Property','Value',...) creates a new FIGBUDDY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before figBuddy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to figBuddy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figBuddy

% Last Modified by GUIDE v2.5 20-Jun-2018 10:20:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figBuddy_OpeningFcn, ...
                   'gui_OutputFcn',  @figBuddy_OutputFcn, ...
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


% --- Executes just before figBuddy is made visible.
function figBuddy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figBuddy (see VARARGIN)

% Choose default command line output for figBuddy
set(hObject,'WindowStyle','normal');

handles.directory='C:\Users\adam.neumann\Documents\GitHub\pc_analysis\figure buddy\presets\';
fn=dir([handles.directory '*.mat']);
fn=struct2cell(fn);
list=fn(1,:);
handles.list=strrep(list,'.mat','');
handles.fn=arrayfun(@(x) [fn{2,x} '\' fn{1,x}],1:size(fn,2),'uniformoutput',false);

handles.hfig=varargin{1};
set(handles.hfig,'WindowStyle','normal');
handles.listbox.String=handles.list;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figBuddy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = figBuddy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox
handles.save_name.String=handles.list{hObject.Value};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function save_name_Callback(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_name as text
%        str2double(get(hObject,'String')) returns contents of save_name as a double


% --- Executes during object creation, after setting all properties.
function save_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=handles.listbox.Value;
load(handles.fn{idx});
set(handles.hfig,'Position',pos);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos=get(handles.hfig,'Position');
save([handles.directory handles.save_name.String],'pos');
new=cellfun(@(x) strcmp(handles.save_name.String,x), handles.list);
if ~any(new)
    handles.list{end+1}=handles.save_name.String;
    handles.fn{end+1}=[handles.directory handles.save_name.String '.mat'];
    handles.listbox.String=handles.list;
end
guidata(hObject, handles);


function update_presets(handles)

