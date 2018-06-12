function varargout = statLTSP_gui(varargin)
% STATLTSP_GUI MATLAB code for statLTSP_gui.fig
%      STATLTSP_GUI, by itself, creates a new STATLTSP_GUI or raises the existing
%      singleton*.
%
%      H = STATLTSP_GUI returns the handle to a new STATLTSP_GUI or the handle to
%      the existing singleton*.
%
%      STATLTSP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STATLTSP_GUI.M with the given input arguments.
%
%      STATLTSP_GUI('Property','Value',...) creates a new STATLTSP_GUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before statLTSP_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to statLTSP_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help statLTSP_gui

% Last Modified by GUIDE v2.5 24-Sep-2017 20:52:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @statLTSP_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @statLTSP_gui_OutputFcn, ...
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



% --- Executes just before statLTSP_gui is made visible.
function statLTSP_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to statLTSP_gui (see VARARGIN)

% Choose default command line output for statLTSP_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);





% UIWAIT makes statLTSP_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = statLTSP_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% First plot: schematic
function plotSchematic(handles)
        loadColors

        %clear(handles.plot_schematic);
        %set(handles.plot_schematic, 'clear');
        b = str2double(get(handles.bound_box, 'String'));
        Prel = max(min(str2double(get(handles.prel_box, 'String')), 0.999), 0.001);
        q = str2double(get(handles.q_box, 'String'));
        
        x = 0:0.005:max(b,1)+0.05;
        y = normpdf(x, Prel*q, sqrt(q^2*Prel*(1-Prel)));
        
        h_rel = plot(handles.plot_schematic, x, y./max(y), '-k', 'linewidth', 3, 'Color', red1);
        hold(handles.plot_schematic, 'on');
        
        y = normpdf(x, b, 0.001);
        h_b = plot(handles.plot_schematic, x, y./max(y), '-k', 'linewidth', 3, 'Color', green2);
        hold(handles.plot_schematic, 'off');
        box(handles.plot_schematic, 'off');
        xlim(handles.plot_schematic, [0 max(b,1)+0.05]);
        ylim(handles.plot_schematic, [0 1.05]);
        
        updateKL(handles);
        updateFlow(handles);
        
        xlabel(handles.plot_schematic, 'postsynaptic potential (mV)', 'fontsize', 10);
        ylabel(handles.plot_schematic, 'norm. freq.', 'fontsize', 10);
        set(handles.plot_schematic, 'TickDir', 'out');
        

function updateKL(handles)        
        KL = @(p,q,t) log(sqrt((q.^2).*p.*(1-p))) + ((t-p.*q).^2)./(2.*((q.^2).*p.*(1-p)));
        %minKL = KL(0.999,1,1);
        minKL = 0;
        klr = KL(str2double(get(handles.prel_box, 'String')), str2double(get(handles.q_box, 'String')), str2double(get(handles.bound_box, 'String'))) + abs(minKL);
        set(handles.kl_box, 'String', num2str((max(klr,0)), 2));
        
function updateFlow(handles)        
        KL = @(p,q,t) log(sqrt((q.^2).*p.*(1-p))) + ((t-p.*q).^2)./(2.*((q.^2).*p.*(1-p)));
        
        cla(handles.plot_flow)
        loadColors
        %hold(handles.plot_flow, 'off');
    
        minp = 0.01;
        minq = 0.01;
        maxp = 1;
        b = str2double(get(handles.bound_box, 'String'));
        maxq = max(b,1)+0.05;
        


        npartsq = 1;
        npartsp = 1;
        incp = 0.01;
        incq = 0.01;

        nq = maxq/npartsq;
        np = maxp/npartsp;
        
        xlim(handles.plot_flow, [0 b+0.05]);

        for u=1:npartsq

            for j=1:npartsp

                Ps = [minp + np*(j-1): incp : maxp - np*(npartsp-j)];
                Qs = [minq + nq*(u-1): incq : maxq - nq*(npartsq-u)];

                KL_M = zeros(length(Ps), length(Qs));
                for z=1:length(Ps)
                    for y=1:length(Qs)
                        KL_M(z,y) = -KL(Ps(z), Qs(y), b);
                    end
                end

                Y = repmat(Ps, size(Qs,2),1)';
                X = repmat(Qs, size(Ps,2),1);

                [U,V] = gradient(KL_M,0.1);

                hold on
                hsl = streamslice(handles.plot_flow, X,Y,U,V, 0.5, 'arrows');
                axis tight
                %hq.MaxHeadSize = 1;
            end
        end

    for j=1:size(hsl,1)
        hsl(j).Color = [0.1 0.1 0.1];
        hsl(j).LineWidth = 0.75;
    end
    
    hold(handles.plot_flow, 'on');    
    if(str2double(get(handles.bound_box, 'String'))>0)
        scatter(handles.plot_flow, str2double(get(handles.bound_box, 'String')), 0.999, 50, green2,'x','Linewidth', 3);
        scatter(handles.plot_flow, str2double(get(handles.q_box, 'String')), str2double(get(handles.prel_box, 'String')), 50, red1,'o','Linewidth', 3);
    else
        plot([0 0], [0 1], 'Color', green2, 'Linewidth', 3);
        plot([0 1], [0 0], 'Color', green2, 'Linewidth', 3);
        scatter(handles.plot_flow, str2double(get(handles.q_box, 'String')), str2double(get(handles.prel_box, 'String')), 50, red1,'o','Linewidth', 3);
    end
    
    
    xlabel(handles.plot_flow, 'q, quantal amp. (mV)', 'Fontsize', 10)
    ylabel(handles.plot_flow, 'P_{rel}, release prob.', 'Fontsize', 10)
    ylim(handles.plot_flow, [0, 1.04]);
    set(handles.plot_flow, 'TickDir', 'out');
    
        


% --- Executes during object creation, after setting all properties.
function prel_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prel_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prel_box_Callback(hObject, eventdata, handles)
% hObject    handle to prel_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prel_box as text
%        str2double(get(hObject,'String')) returns contents of prel_box as a double
handles.data.prel_v = max(min(str2double(get(hObject, 'String')),0.999),0.001);
if isnan(handles.data.prel_v)
    errordlg('Input must be a number','Error');
end
set(handles.prel_box, 'String', num2str(handles.data.prel_v));
set(handles.prel_slider, 'Value', handles.data.prel_v);

plotSchematic(handles);

% Save the new prel_box value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function q_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to q_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function q_box_Callback(hObject, eventdata, handles)
% hObject    handle to q_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of q_box as text
%        str2double(get(hObject,'String')) returns contents of q_box as a double
handles.data.q_v = str2double(get(hObject, 'String'));
set(handles.q_slider, 'Value', handles.data.q_v);
if isnan(handles.data.q_v)
    errordlg('Input must be a number','Error');
end

plotSchematic(handles);
guidata(hObject,handles)



% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

set(handles.prel_box, 'String', 0.5);
set(handles.q_box, 'String', 0.7);
set(handles.bound_box, 'String', 1);
set(handles.kl_box, 'String', 0);

set(gcf, 'Name', 'Statistical Long-term Synaptic Plasticity (statLTSP)');
set(gcf, 'Resize', 'off');

handles.data.prel_v = str2double(get(handles.prel_box, 'String'));
handles.data.q_v = str2double(get(handles.q_box, 'String'));
handles.data.bound_v = str2double(get(handles.bound_box, 'String'));

plotSchematic(handles);

% Update handles structure
guidata(handles.figure1, handles);



function bound_box_Callback(hObject, eventdata, handles)
% hObject    handle to bound_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.bound_v = min(str2double(get(hObject,'String')),3);
set(handles.bound_box, 'String', num2str(handles.data.bound_v));
set(handles.q_slider, 'Max', max(handles.data.bound_v,1));
set(handles.bound_slider, 'Value', handles.data.bound_v);
plotSchematic(handles);


% --- Executes during object creation, after setting all properties.
function bound_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bound_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function prel_slider_Callback(hObject, eventdata, handles)
% hObject    handle to prel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.data.prel_v = max(min(get(hObject,'Value'),0.999),0.001);
set(handles.prel_box, 'String', num2str(handles.data.prel_v));
plotSchematic(handles);


% --- Executes during object creation, after setting all properties.
function prel_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function q_slider_Callback(hObject, eventdata, handles)
% hObject    handle to q_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.q_v = get(hObject,'Value');
set(handles.q_box, 'String', num2str(handles.data.q_v));
plotSchematic(handles);


% --- Executes during object creation, after setting all properties.
function q_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to q_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function bound_slider_Callback(hObject, eventdata, handles)
% hObject    handle to bound_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.bound_v = get(hObject,'Value');
set(handles.bound_box, 'String', num2str(handles.data.bound_v));
set(handles.kl_box, 'String', num2str(handles.data.bound_v));
set(handles.q_slider, 'Max', max((handles.data.bound_v),1));

plotSchematic(handles);
updateKL(handles);


% --- Executes during object creation, after setting all properties.
function bound_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bound_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
