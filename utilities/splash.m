function varargout = splash(varargin)
% SPLASH MATLAB code for splash.fig
%      SPLASH, by itself, creates a new SPLASH or raises the existing
%      singleton*.
%
%      H = SPLASH returns the handle to a new SPLASH or the handle to
%      the existing singleton*.
%
%      SPLASH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLASH.M with the given input arguments.
%
%      SPLASH('Property','Value',...) creates a new SPLASH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before splash_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to splash_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help splash

% Last Modified by GUIDE v2.5 28-Apr-2015 18:07:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @splash_OpeningFcn, ...
                   'gui_OutputFcn',  @splash_OutputFcn, ...
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


% --- Executes just before splash is made visible.
function splash_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to splash (see VARARGIN)

%pixels
set( handles.figure1, ...
    'Units', 'pixels' );

%get your display size
screenSize = get(0, 'ScreenSize');

position = get( handles.figure1, ...
    'Position' );

position(1) = (screenSize(3)-position(3))/2;
position(2) = (screenSize(4)-position(4))/2;

set( handles.figure1, ...
    'Position', position );

% Choose default command line output for splash
handles.output = hObject;

try
    img = imread( './utilities/logo.jpg' );
catch err
    img = imread( 'logo.jpg' );
end

imshow( img );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes splash wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = splash_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in closeButton.
function closeButton_Callback(hObject, eventdata, handles)
% hObject    handle to closeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;
%icaoberg 9/28/2013
%hide toolbar from this release
%toolbar;


% --------------------------------------------------------------------
function web_Callback(hObject, eventdata, handles)
% hObject    handle to web (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
