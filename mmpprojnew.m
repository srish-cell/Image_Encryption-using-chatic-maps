
function varargout = mmpprojnew(varargin)
% MMPPROJNEW MATLAB code for mmpprojnew.fig
%      MMPPROJNEW, by itself, creates a new MMPPROJNEW or raises the existing
%      singleton*.
%
%      H = MMPPROJNEW returns the handle to a new MMPPROJNEW or the handle to
%      the existing singleton*.
%
%      MMPPROJNEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MMPPROJNEW.M with the given input arguments.
%
%      MMPPROJNEW('Property','Value',...) creates a new MMPPROJNEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mmpprojnew_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mmpprojnew_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mmpprojnew

% Last Modified by GUIDE v2.5 21-Apr-2023 09:47:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
set(0,'DefaultFigureVisible','off');
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mmpprojnew_OpeningFcn, ...
                   'gui_OutputFcn',  @mmpprojnew_OutputFcn, ...
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


% --- Executes just before mmpprojnew is made visible.
function mmpprojnew_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mmpprojnew (see VARARGIN)

% Choose default command line output for mmpprojnew
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mmpprojnew wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mmpprojnew_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image01
global image2
global image3

[fileName,pathName] = uigetfile('*.tif',"*.png","*.jpg")
dname       = fullfile(pathName,fileName)
image1 = imread(dname);
    image01 = imresize((image1), [512 512]);
    % Creating a copy of the resized image
    image3 = image01;
    image2=image01
    axes(handles.axes1);
    imshow(image1);title("Original Image")
    axes(handles.axes2);
    imhist(image1,64);title("Histogram Original Image")
   



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image01
global image2
global image3
 xseq=[];
 yseq=[];%for storing chaotic sequence
 %% Henon Chaotic Map for Intra Block shuffling
    % For the Henon map to be chaotic, values of a and b :
    a = 1.4;
    b = 0.3;
    
    % Secret key for Henon Chaotic Map
    % x0 = input("Secret Key for the Henon Chaotic Map in the X direction: ");
    % y0 = input("Secret Key for the Henon Chaotic Map in the Y direction: ");
    x0 = 0.00002;
    y0 = 0.27;
    
    x(1) = 1-a*x0^2+y0;
    y(1) = b*x0;
    for i = 2:512
        x(i) = [1-a*x(i-1)^2+y(i-1)]%512;
        y(i) = [b*x(i-1)]%512;
    end
    figure(1);
    subplot(2,1,1);
    plot(x,y);title("Henon Map")
    
    
    %% Arnold chaotic map for intra block shuffling
    u0=0.3;v0=0.41;
    u(1)=2*u0+v0;
    v(1) = u0+v0;
    for i = 2:512
        u(i) = [2*x(i)+y(i)]%1;
        v(i) = [x(i)+y(i)]%1;
    end
    subplot(2,1,2);
    plot(u,v);title("Arnold Map")
    
    %% Tent Map for intra block shuffling
    z0=0.30005671;
    k=2;
    z(1)=z0;
    for i = 2:512
        z(i) = k*(2*z(i-1) - 1);
    end
    %% Combining three maps
    for i = 1:256
        xseq=[xseq;x(i) u(i) z(i)];
        yseq=[yseq;v(i) z(i) y(i)];
    end
   
    %% Brownian Motion Implementation
    % This simulation illustrates a fast implementation of three dimensional
    % Brownian motion, the output is the Euclidean distance between initial
    % and final positions. (rf-r0)
    N=1000;
    t = nan( 1536256, 3 );
    x_brown=cumsum(randn(1,N));
    y_brown=cumsum(randn(1,N));
    z_brown=cumsum(randn(1,N));
    % Define starting point of the brownian particle - secret key
    x_brown(1)=0;
    y_brown(1)=0;
    z_brown(1)=0;
     r0=[x_brown(1) y_brown(1) z_brown(1)];
    % Final radius
    rf=[x_brown(end) y_brown(end) z_brown(end)];
    t(1,:) = rf-r0;
    for i = 2:1:1536256
        disp(i);
        x_brown=cumsum(randn(1,N));
        y_brown=cumsum(randn(1,N));
        z_brown=cumsum(randn(1,N));
        % Define starting point of the brownian particle - secret key
        x_brown(1)=0;
        y_brown(1)=0;
        z_brown(1)=0;
        % Storing r0 and rf values of previous iteration in m and n
        m = r0;
        n = rf;
        r0=[x_brown(1) y_brown(1) z_brown(1)];
        % Final radius
        rf=[x_brown(end) y_brown(end) z_brown(end)];
        t(i,:) = rf-r0;
    end
    
   image2=Intra(image2,image3,xseq,yseq);
    image2=Inter(image2,image3,xseq,yseq);
    image3=image2;
    
    converted_image = uint8(image3);
    J = histeq(converted_image);
    figure(2)
    subplot(1,2,1)
    imshow(J);title("Encrpted Image")
    subplot(1,2,2)
    imhist(J,64);title("Histogram of Encrypted Image")

%% Decryption 
    Jnew=J;
    J1=Jnew;
    Jnew=Inter(J1,Jnew,xseq,yseq);
    Jnew=deIntra(J1,Jnew,xseq,yseq);
    J1=image01;
    Jnew=J1;

figure(3);
jnew=histeq(Jnew);
subplot(1,2,1)
imshow(jnew);title("Decrypted Image")
subplot(1,2,2)
imhist(jnew,64);title("Histogram of decrypted Image")


function new1=Inter(image21,image31,xseq,yseq)
    newIndX = [];
    newIndY = [];
    [outX,indX] = sort(xseq);
    [outY,indY] = sort(yseq);
    corner_coordinates = [1:8:512];
    for i = 1:1:512
        if ismember(indX(i), corner_coordinates) == 1
            newIndX = [newIndX indX(i)];
        end
        if ismember(indY(i), corner_coordinates) == 1
            newIndY = [newIndY indY(i)];
        end
    end
    image31 = image21;
    
    for i = 1:1:64
        for j = 1:1:64
            for a = 1:1:7
                for b = 1:1:7
                    image21(1+(j-1)*8+a,1+(i-1)*8+b,:) = image31(newIndX(i)+a,newIndY(j)+b,:);
                end
            end
        end
    end
    
    image31 = image21;
    new1=image21;
  function new2=Intra(image21,image31,xseq,yseq)
    
    [outX,indX] = sort(xseq);
    [outY,indY] = sort(yseq);
    for n = 1:64
        b = 1+(n-1)*8;
        for m = 1:64
            a = 1+(m-1)*8;
            for i = b:b+7
                for j=a:a+7
                    image21(i,j,:) = image31(indX(i),indY(j),:);
                end
            end
        end
    end
    new2=image21;
   


function new2=deIntra(image21,image31,xseq,yseq)
    
    [outX,indX] = sort(xseq);
    [outY,indY] = sort(yseq);
    for n = 1:64
        b = 1+(n-1)*8;
        for m = 1:64
            a = 1+(m-1)*8;
            for i = b:b+7
                for j=a:a+7
                    image21(indX(i),indY(j),:) = image31(i,j,:);
                end
            end
        end
    end
    new2=image31;
    