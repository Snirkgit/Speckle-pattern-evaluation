function varargout = RealDICSimGui(varargin)
% REALDICSIMGUI MATLAB code for RealDICSimGui.fig
%      REALDICSIMGUI, by itself, creates a new REALDICSIMGUI or raises the existing
%      singleton*.
%
%      H = REALDICSIMGUI returns the handle to a new REALDICSIMGUI or the handle to
%      the existing singleton*.
%
%      REALDICSIMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REALDICSIMGUI.M with the given input arguments.
%
%      REALDICSIMGUI('Property','Value',...) creates a new REALDICSIMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RealDICSimGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RealDICSimGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RealDICSimGui

% Last Modified by GUIDE v2.5 31-May-2018 21:13:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RealDICSimGui_OpeningFcn, ...
                   'gui_OutputFcn',  @RealDICSimGui_OutputFcn, ...
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


% --- Executes just before RealDICSimGui is made visible.
function RealDICSimGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RealDICSimGui (see VARARGIN)
global Completed
Completed = 0 ;
handles.CData.Data = [0 0;0 0] ;
handles.SbData.Data = [0 0;0 0;0 0;0 0] ;

% Choose default command line output for RealDICSimGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RealDICSimGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RealDICSimGui_OutputFcn(hObject, eventdata, handles) 
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
global SubData ; global SubXSiz ; global SubYSiz ;
global OrigIm ; global WarpIm ; global RgbIm ; global RGmap ;
global AnaField ; global RecField ; global ErrField ; global Completed ;
handles.Load.ForegroundColor = [0 0 0] ;
% T is the Affine Deformation Matrix
% Principal Diagonal - Axial Tenstion
% Secondary Diagonal - Shear deformaion
% for more information about T, search for affine2d() function.
T = [(100 + str2num(handles.XTen.String))/100 str2num(handles.Shr.String)/100 0; str2num(handles.Shr.String)/100 (100 + str2num(handles.YTen.String))/100 0; 0 0 1] ;

% OrigIm is the input image - as grayscale format, that you wish to check it's speckle pattern
%OrigIm = imread('Fibers Pattern.tif') ;

% tform is a image transformation Object, it is to be used in the imwarp() function.
tform = affine2d(T) ;

% WarpIm is the warped image as grayscale 2D matrix
WarpIm = imwarp(OrigIm,tform) ;

% The next lines include the most important input parameters
% SubXSiz - the size of each rectengular Subset set to be correlated in x direction
% SubYSiz - the size of each rectengular Subset set to be correlated in y direction
% X_Spacing - the size of each step to be taken from the left end of the Subset in x direction
% Y_Spacing - the size of each step to be taken from the upper end of the Subset in y direction
% Extend - this parameter allow you to extend the searching areas to be correlated.
% The value of this parameter indicates the extension of the searching zone in each direction - (x & y)
% RGThreshold - The threshold value of the absolute error for the Red Green image
SubXSiz = str2num(handles.XSbSz.String) ;  SubYSiz = str2num(handles.YSbSz.String) ;
X_Spacing = str2num(handles.XSpSz.String) ;    Y_Spacing = str2num(handles.YSpSz.String) ;
Extend = str2num(handles.Extnd.String)  ;        RGThreshold = 2 ;

% The next lines include all the parameters that need to be intialize, and not to be changed
[OrigImYSiz , OrigImXSiz] = size(OrigIm) ;
AvError = [0 0];    SubsetAbsError = [0 0];
RgbIm = cat(3, OrigIm, OrigIm, OrigIm) ;
X_Pix = 1; Y_Pix = 1;       Runs = 0;  SubData = [];

% The main loops to choose subsets and correlate them are starting here
% The first while loop runs over y direction and the second runs over x direction
% X_Pix and Y_Pix indicates the current subset's top left corner pixels
% U_Former and V_Former are values that helps the algorithm to predict the
% next location to correlate, the appearance of them now is just to give them the initialize value, and not to be changed
while Y_Pix + SubYSiz - 1 <= OrigImYSiz
    U_Former = 0 ; V_Former = 0 ;
    while 1
        
        i0 = Y_Pix ; i1 = Y_Pix - 1 + SubYSiz ;
        j0 = X_Pix ; j1 = X_Pix - 1 + SubXSiz ;
        
        % MidSubY and MidSubX are the middle pixels of the current subset
        % u is the analytical displacement of the middle of the current subset according to the analytical displacement field defined by T
        MidSubY = floor((i1+i0)/2) ;
        MidSubX = floor((j0+j1)/2) ;
        R = [MidSubX ; MidSubY ; 0] ;
        Rnew = T*[MidSubX ; MidSubY ; 0] ;
        u = Rnew - R ;
        
        % Cut the current pair of subsets to be correlated according to the prediction
        OrigImSub = OrigIm( i0 : i1 , j0 : j1 ) ;
        WarpImSub = WarpIm( i0 + V_Former : i1 + V_Former , j0 + U_Former : j1 + U_Former ) ;
        
        % The Extend feature described above is applied in the next lines
        % Extend the Subset for wider Search window
        m0 = OrigIm ;
        m1 = WarpIm ;
        
        [M , N]=size(m0);
                        
        dm0 = OrigImSub ;
        dm1 = WarpImSub ;
                
        % In this following loop, we stated the conditions for extending the image
        while 1
            Extend_y_0 = Extend ; Extend_y_1 = Extend ;
            Extend_x_0 = Extend ; Extend_x_1 = Extend ;
            
            if Extend > i0 ,          Extend_y_0 = i0 ;         end
            if Extend + i1 > M ,    Extend_y_1 = M - i1 ;	end
            if Extend  > j0 ,         Extend_x_0 = j0 ;        end
            if Extend + j1 > N ,    Extend_x_1 = N - j1 ;   end
            
            % Condition 1 - If the user chose Extend > 0, the subset will be extended by Extend
            dm1Extended_Size = m1(1 + i0 - Extend_y_0 : i1 + Extend_y_1 , 1 + j0 - Extend_x_0 : j1 + Extend_x_1 ) ;
            [M1 , N1] = size(dm1) ;
            dm1Extended = zeros(M1 + 2*Extend , N1 + 2*Extend  ,'uint8') ;
            dm1Extended(2 + Extend - Extend_y_0 :  M1 + 2*Extend - (Extend - Extend_y_1)  , 2 + Extend - Extend_x_0 : ...
                N1 + 2*Extend - (Extend - Extend_x_1)) = dm1Extended_Size;
            
            % Add zeros to dm1 just to make sure both of the subsets having the same size, and thus can be correlated by decorr2a() function
            dm0Extended = zeros(M1 + 2*Extend , N1 + 2*Extend ,'uint8') ;
            dm0Extended(1 + Extend : M1 + Extend, 1 + Extend : N1 + Extend) = dm0 ;
            
            % Codition 2 -  If even one of the subsets have not enough information, i.e. they are all zeros - extend the subset by 10
            if any(any(dm0Extended)) * any(any(dm1Extended))
                    if handles.FFTCC.Value == 1
                    % Option 1 - Cross-Correlation by Franck
                    [du,dv,dc] = decorr2a(dm0Extended-mean2(dm0Extended),dm1Extended-mean2(dm1Extended)) ;
                    end
                    if handles.NCC.Value == 1
                        if all(all(dm0 == dm0(1,1)))
                            dm0(1,1) = dm0(1,1) - 1 ;
                        else
                            % Option 2 - Cross-Correlation. Better results, but slower for big Subsets. Added in 29/5
                            c = normxcorr2(dm0,dm1Extended) ;
                            [ypeak, xpeak] = find(c==max(c(:))) ;
                            du = xpeak - (N1 + Extend) ;
                            dv = ypeak - (M1 + Extend) ;
                        end
                    end
                    if handles.FFTNCC.Value == 1
                            % Option 3 - template_matching function by FFT NCC
                            c = template_matching(dm0,dm1Extended) ;
                            [ypeak, xpeak] = find(c==max(c(:))) ;
                            du = xpeak  - (N1/2 + Extend) ;
                            dv = ypeak  - (M1/2 + Extend) ;
                    end
                break
            else
                Extend = Extend + 10 ;
            end
        end
        
        du = du(1) ;
        dv = dv(1) ;
        
        % This variable contains the correlated displacement
        Rcorrelated = [du ; dv ; 0]  ;
        
        % Compute the errors of the correlation, this is not the final results.
        % The results of AvError will be given after the main while loop.
        
        if ~isnan(du) && ~isnan(dv)
            
            %{
            % The AvError can be the either the Absolute Average Error or the Average Relative Error.
            % The default output error is the Absolute, but one can choose to use the Relative if he wish to do so.
            AvError(1) = AvError(1)  + u(1) - Rcorrelated(1) ; % AvError(1) = AvError(1) +  100*((u(1) - Rcorrelated(1))/(u(1)) );
            AvError(2) = AvError(2) + u(2) - Rcorrelated(2) ; %AvError(2) = AvError(2) + 100*((u(2) - Rcorrelated(2))/(u(2)) );
            %}
            
            % This error is the absolute error - # of deviated pixels.
            SubsetAbsError(1) = u(1) - Rcorrelated(1) ;
            SubsetAbsError(2) = u(2) - Rcorrelated(2) ;
            
            
            % The following if statements, determining if the subset will be colored by Red of Green according to the RGThreshold from the user.
            if norm(SubsetAbsError) < RGThreshold
                RgbIm(i0:i1,j0:j1,2) = RgbIm(i0:i1,j0:j1,2) + 30/(floor(SubXSiz/X_Spacing)+1) ;
            end
            if norm(SubsetAbsError) >= RGThreshold
                RgbIm(i0:i1,j0:j1,1) = RgbIm(i0:i1,j0:j1,1) + 30/(floor(SubXSiz/X_Spacing)+1);
            end
            
            % The variable SubData include all the data from the simulation.
            % It contains the subset position, his analytical displacement and  his correlated displacement.
            SubData(Runs+1,1:6) = [MidSubX ; MidSubY ; u(1) ; u(2) ; du ; dv] ;
            
        end
        Runs = Runs + 1 ;
        handles.Load.Visible = 'on' ;
        if mod(Runs,500) == 0
            handles.Load.String = ['Computing  ' num2str(100*Runs/((M/X_Spacing)*(N/Y_Spacing))) '%'];
            pause(0.0000001)
        end
        
        %{ 
        % Using the former steps to predict the next one.
        du = du + U_Former ; dv = dv + V_Former ;
        U_Former = floor(du) ; V_Former = -floor(dv) ;
        %}        
        
        % Updating the parameters before the next run
        X_Pix = X_Pix + X_Spacing ;
        if X_Pix + SubXSiz - 1 > OrigImXSiz
            X_Pix = 1 ; Y_Pix = Y_Pix + Y_Spacing ;
            break
        end
        
    end
      
end

handles.Load.String = ['Computing   100%'];
handles.Load.ForegroundColor = [0 0.7 0.5] ;

disp(' ')
disp('-------Correlation Data-------')
disp(' ')

% Final Error calculation & display
% AvError = AvError/Runs ;
SubData(:,7) = SubData(:,3) - SubData(:,5) ;
SubData(:,8) = SubData(:,4) - SubData(:,6) ;

AvError = [mean(SubData(:,7)) mean(SubData(:,8))] ;
NormError = norm(AvError) ;
disp('The Average Error is:')
disp(AvError)
%disp(NormError)

% Calculate Standard Deviation
Std = [std(SubData(:,7)) std(SubData(:,8))] ;
disp('The Standard Deviation is:')
disp(Std)

CorData = [AvError(1) AvError(2) ; Std(1) Std(2)] ;
handles.CData.Data = CorData ;
handles.text29.Visible = 'on' ;
handles.CData.Visible = 'on' ;

% Showing all the visual results
%{
OriFig = figure('name','Original Image','NumberTitle','off') ;
imshow(OrigIm)
title('Original Image')
%}

%DefFig = figure('name', 'Deformed Image','NumberTitle','off') ;
axes(handles.axes2)
imshow(WarpIm) ; % hold on added in 29/5
title('Deformed Image')

%RGFig = figure('name', 'Red-Green Error Map','NumberTitle','off') ;
axes(handles.axes1)
imshow(OrigIm) ; hold on% hold on added in 29/5
title('Original Image')

RGmap = imshow(RgbIm) ;
if handles.checkbox1.Value == 1
    RGmap.Visible = 'on' ;
else
    RGmap.Visible = 'off' ;
end

AnaField = quiver(SubData(:,1),SubData(:,2),SubData(:,3),SubData(:,4),0,'Color','w') ;
if handles.checkbox2.Value == 1
    AnaField.Visible = 'on' ;
else
    AnaField.Visible = 'off' ;
end

RecField = quiver(SubData(:,1),SubData(:,2),SubData(:,5),SubData(:,6),0,'Color','[1 0.6 0.9]') ;
if handles.checkbox3.Value == 1
    RecField.Visible = 'on' ;
else
    RecField.Visible = 'off' ;
end

ErrField = quiver(SubData(:,1),SubData(:,2),SubData(:,7),SubData(:,8),0,'Color','[0 1 1]') ;
if handles.checkbox4.Value == 1
    ErrField.Visible = 'on' ;
else
    ErrField.Visible = 'off' ;
end

handles.PresDat.Visible = 'on' ;

Completed = 1 ;
%%% - End - %%%


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SubData ; global SubXSiz ; global SubYSiz ;
global WarpIm ; global RgbIm ; global OrigIm ; global Completed ;
global RGmap ; global AnaField ; global RecField ;  global ErrField ; 

if Completed == 1
    if strcmp(handles.figure1.SelectionType,'alt') == 1
        [XLoc,YLoc] = ginput(1) ;
        SubLoc = SubData(dsearchn(SubData(:,1:2),[XLoc YLoc]),1:8) ;
        disp('-------Subset Data-------')
        disp(' ')
        disp('Subset Position:')
        disp(SubLoc(1:2))
        disp('Analytical Displacement:')
        disp(SubLoc(3:4))
        disp('Recovery Displacement:')
        disp(SubLoc(5:6))
                
        SubDat = [SubLoc(1:2) ; SubLoc(3:4) ; SubLoc(5:6) ; SubLoc(7:8)] ;
        handles.SbData.Data = SubDat ;
        handles.SbData.Visible = 'on' ;
        handles.text30.Visible = 'on' ;
        
        %ginput(0) ;
        axes(handles.axes1)
        rectangle('Position',[SubLoc(1)-SubYSiz/2,SubLoc(2)-SubXSiz/2,SubYSiz,SubXSiz],...
                'Curvature',[0.1,0.1],...
                'LineWidth',2,'LineStyle','-.','EdgeColor','yellow') ;
        
        axes(handles.axes2)
        rectangle('Position',[SubLoc(1)-SubYSiz/2+SubLoc(3),SubLoc(2)+SubLoc(4)-SubXSiz/2,SubYSiz,SubXSiz],...
                'Curvature',[0.1,0.1],...
                'LineWidth',2,'LineStyle','-.','EdgeColor','white') ;
        rectangle('Position',[SubLoc(1)-SubYSiz/2+SubLoc(5),SubLoc(2)+SubLoc(6)-SubXSiz/2,SubYSiz,SubXSiz],...
                'Curvature',[0.1,0.1],...
                'LineWidth',2,'LineStyle','-.','EdgeColor','[1 0.411 0.705]') ;
        
    end

    if strcmp(handles.figure1.SelectionType,'extend') == 1
        
        axes(handles.axes2)
        imshow(WarpIm) ; % hold on added in 29/5
        title('Deformed Image')

        axes(handles.axes1)
        imshow(OrigIm) ; hold on % hold on added in 29/5
        title('Original Image')

        RGmap = imshow(RgbIm) ;
        if handles.checkbox1.Value == 1
            RGmap.Visible = 'on' ;
        else
            RGmap.Visible = 'off' ;
        end

        AnaField = quiver(SubData(:,1),SubData(:,2),SubData(:,3),SubData(:,4),0,'Color','w') ;
        if handles.checkbox2.Value == 1
            AnaField.Visible = 'on' ;
        else
            AnaField.Visible = 'off' ;
        end

        RecField = quiver(SubData(:,1),SubData(:,2),SubData(:,5),SubData(:,6),0,'Color','[1 0.6 0.9]') ;
        if handles.checkbox3.Value == 1
            RecField.Visible = 'on' ;
        else
            RecField.Visible = 'off' ;
        end

        ErrField = quiver(SubData(:,1),SubData(:,2),SubData(:,7),SubData(:,8),0,'Color','[0 1 1]') ;
        if handles.checkbox4.Value == 1
            ErrField.Visible = 'on' ;
        else
            ErrField.Visible = 'off' ;
        end

        handles.PresDat.Visible = 'on' ;
    end
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RGmap ; global Completed
if Completed == 1
    if handles.checkbox1.Value == 1
        RGmap.Visible = 'on' ;
    else
        RGmap.Visible = 'off' ;
    end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global AnaField ; global Completed
if Completed == 1
    if handles.checkbox2.Value == 1
        AnaField.Visible = 'on' ;
    else
        AnaField.Visible = 'off' ;
    end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Completed ; global RecField
if Completed == 1
    if handles.checkbox3.Value == 1
        RecField.Visible = 'on' ;
    else
        RecField.Visible = 'off' ;
    end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Completed ; global ErrField
if Completed == 1
    if handles.checkbox4.Value == 1
        ErrField.Visible = 'on' ;
    else
        ErrField.Visible = 'off' ;
    end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox4



function XSbSz_Callback(hObject, eventdata, handles)
% hObject    handle to XSbSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XSbSz as text
%        str2double(get(hObject,'String')) returns contents of XSbSz as a double


% --- Executes during object creation, after setting all properties.
function XSbSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XSbSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YSbSz_Callback(hObject, eventdata, handles)
% hObject    handle to YSbSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YSbSz as text
%        str2double(get(hObject,'String')) returns contents of YSbSz as a double


% --- Executes during object creation, after setting all properties.
function YSbSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YSbSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XTen_Callback(hObject, eventdata, handles)
% hObject    handle to XTen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OrigIm ; global Browsed ;
if Browsed == 1 ;
    cla(handles.axes1,'reset')
    axes(handles.axes1)
    imshow(OrigIm)
    title('Original Image')

    T = [(100 + str2num(handles.XTen.String))/100 str2num(handles.Shr.String)/100 0; str2num(handles.Shr.String)/100 (100 + str2num(handles.YTen.String))/100 0; 0 0 1] ;
    tform = affine2d(T) ;
    WarpIm = imwarp(OrigIm,tform) ;
    cla(handles.axes2,'reset')
    axes(handles.axes2)
    imshow(WarpIm)
    title('Deformed Image')
end
% Hints: get(hObject,'String') returns contents of XTen as text
%        str2double(get(hObject,'String')) returns contents of XTen as a double


% --- Executes during object creation, after setting all properties.
function XTen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XTen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YTen_Callback(hObject, eventdata, handles)
% hObject    handle to YTen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OrigIm ; global Browsed ;
if Browsed == 1 ;
    cla(handles.axes1,'reset')
    axes(handles.axes1)
    imshow(OrigIm)
    title('Original Image')

    T = [(100 + str2num(handles.XTen.String))/100 str2num(handles.Shr.String)/100 0; str2num(handles.Shr.String)/100 (100 + str2num(handles.YTen.String))/100 0; 0 0 1] ;
    tform = affine2d(T) ;
    WarpIm = imwarp(OrigIm,tform) ;
    cla(handles.axes2,'reset')
    axes(handles.axes2)
    imshow(WarpIm)
    title('Deformed Image')
end
% Hints: get(hObject,'String') returns contents of YTen as text
%        str2double(get(hObject,'String')) returns contents of YTen as a double


% --- Executes during object creation, after setting all properties.
function YTen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YTen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shr_Callback(hObject, eventdata, handles)
% hObject    handle to Shr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OrigIm ; global Browsed ;
if Browsed == 1 ;
    cla(handles.axes1,'reset')
    axes(handles.axes1)
    imshow(OrigIm)
    title('Original Image')

    T = [(100 + str2num(handles.XTen.String))/100 str2num(handles.Shr.String)/100 0; str2num(handles.Shr.String)/100 (100 + str2num(handles.YTen.String))/100 0; 0 0 1] ;
    tform = affine2d(T) ;
    WarpIm = imwarp(OrigIm,tform) ;
    cla(handles.axes2,'reset')
    axes(handles.axes2)
    imshow(WarpIm)
    title('Deformed Image')
end
% Hints: get(hObject,'String') returns contents of Shr as text
%        str2double(get(hObject,'String')) returns contents of Shr as a double


% --- Executes during object creation, after setting all properties.
function Shr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Extnd_Callback(hObject, eventdata, handles)
% hObject    handle to Extnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Extnd as text
%        str2double(get(hObject,'String')) returns contents of Extnd as a double


% --- Executes during object creation, after setting all properties.
function Extnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Extnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XSpSz_Callback(hObject, eventdata, handles)
% hObject    handle to XSpSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XSpSz as text
%        str2double(get(hObject,'String')) returns contents of XSpSz as a double


% --- Executes during object creation, after setting all properties.
function XSpSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XSpSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YSpSz_Callback(hObject, eventdata, handles)
% hObject    handle to YSpSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YSpSz as text
%        str2double(get(hObject,'String')) returns contents of YSpSz as a double


% --- Executes during object creation, after setting all properties.
function YSpSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YSpSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key release with focus on figure1 and none of its controls.
function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
global Completed ; global SubData ;
if strcmp(eventdata.Character,'h') == 1
     if Completed == 1;
        x_hist = figure('name', 'X Error Histogram','NumberTitle','off') ;
        figure(x_hist.Number) ; histogram(SubData(:,7)) ; title('X Error Histogram') ;
        %{
        Xlm = get(gca,'XLim') ;         Ylm = get(gca,'YLim') ; 
        Xtxt = {['Average Error: ', num2str(AvError(1)), '  Pixels' ] , ['Standard Deviation: ', num2str(Std(1)), '  Pixels'] } ;
        text(Xlm(2)-(Xlm(2)-Xlm(1))*0.46,Ylm(2)/1.1,Xtxt,'Color',[0.9 0.1 0])
        %}
        xlabel('Subset X Error [Pixels]') ;      ylabel('Number  of  Subsets  [-]') ;
        
        y_hist = figure('name', 'Y Error Histogram','NumberTitle','off') ;
        figure(y_hist.Number) ; histogram(SubData(:,8)) ; title('Y Error Histogram') ;
        %{
        Xlm = get(gca,'XLim') ;         Ylm = get(gca,'YLim') ; 
        Ytxt = {['Average Error: ', num2str(AvError(2)), '  Pixels' ] , ['Standard Deviation: ', num2str(Std(2)), '  Pixels'] } ;
        text(Xlm(2)-(Xlm(2)-Xlm(1))*0.46,Ylm(2)/1.1,Ytxt,'Color',[0.9 0.1 0])
        %}
        xlabel('Subset Y Error [Pixels]') ;      ylabel('Number  of  Subsets  [-]') ;
     end
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OrigIm ; global WarpIm ; global Completed ; global Browsed ;
Browsed = 0 ;
[Filename, Pathname] = uigetfile('*.tif','File Selector') ;
if Filename ~= 0
    name = strcat(Pathname,Filename) ;
    OrigIm = imread(name) ;
    size(size(OrigIm))
    MtchSz = size(size(OrigIm)) == [1 2] ;
    if MtchSz(2) == 0
        OrigIm = rgb2gray(OrigIm(:,:,1:3)) ;
    end
    handles.ImName.String = Filename ;
    cla(handles.axes1,'reset')
    axes(handles.axes1)
    imshow(OrigIm)
    title('Original Image')

    T = [(100 + str2num(handles.XTen.String))/100 str2num(handles.Shr.String)/100 0; str2num(handles.Shr.String)/100 (100 + str2num(handles.YTen.String))/100 0; 0 0 1] ;
    tform = affine2d(T) ;
    WarpIm = imwarp(OrigIm,tform) ;
    cla(handles.axes2,'reset')
    axes(handles.axes2)
    imshow(WarpIm)
    title('Deformed Image')
    
    handles.CData.Data = [0 0;0 0] ;
    handles.SbData.Data = [0 0;0 0;0 0;0 0] ;
    handles.Load.String = 'Computing';
    handles.Load.ForegroundColor = [0 0 0] ;
    Completed = 0 ;
    Browsed = 1 ;
end


% --- Executes on button press in About.
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    msgbox('Developed by Aviya Shmuel & Snir Koska (c)','About','help');


% --- Executes on button press in FFTNCC.
function FFTNCC_Callback(hObject, eventdata, handles)
% hObject    handle to FFTNCC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FFTNCC
