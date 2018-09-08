function adam_gui

global selected_vertex;
%Initialize variables to some bogus values so that they can be set by subfunctions
selected_vertex=-1;
DataFig=-1;
ProjFig=-1;
Proj3DFig=-1;
OCCFig=-1;
OCC=-1;
AOtext='Select AO';
OCCindex=1;
nOCC=0;
Proj3DFigPlot=-1;
axislim=0;
OCCmenu={'OCC1'};
OCCOffset=-1;
current_view=1;
nfac=1;
nvert=1;
TR=0;
popmenu=cell(1,1);
nAO=0;
VN=0;
tlist=0;
vlist=0;
vlistR=0;
FT=0;
angles=0;
im=0;
E=0;
E0=0;
up=0;
Date=0;
PixScale=0;
Filename=0;
RotAngle=0;
min_tim=0;
iniFileName='';
iniPathName='';
Albedo='';
Beta=0;
Lambda=0;
Period=0;

%cd /home/matvii/Storage/Data/final2/7_Iris_4268/;
LoadIni;
%keyboard

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[360,500,450,285]);
hcheckrot = uicontrol('Style','pushbutton','String','Rotate',...
    'Position',[315,260,40,15],...
    'Callback',{@rotbutton_Callback});
hsundir=uicontrol('Style','pushbutton','String','Sun dir','Position',[360,260,40,15], 'Callback',{@hsundir_Callback});
%  Construct the components.
hload = uicontrol('Style','pushbutton','String','Load',...
    'Position',[315,220,40,15],...
    'Callback',{@loadbutton_Callback});
hloadmesh = uicontrol('Style','pushbutton','String','Load mesh',...
    'Position',[360,220,40,15],...
    'Callback',{@loadmeshbutton_Callback});
hloadOCC=uicontrol('Style','pushbutton','String','Load OCC',...
    'Position',[360,190,40,15],...
    'Callback',{@loadOCCbutton_Callback});
hsave = uicontrol('Style','pushbutton','String','Save mesh',...
    'Position',[315,190,40,15],...
    'Callback',{@savebutton_Callback});
hcontour = uicontrol('Style','pushbutton',...
    'String','Projections',...
    'Position',[315,145,40,15],...
    'Callback',{@projbutton_Callback});
hcheckfit = uicontrol('Style','pushbutton',...
    'String','Check fit',...
    'Position',[360,145,40,15],...
    'Callback',{@fitbutton_Callback});
htext = uicontrol('Style','text','String',AOtext,...
    'Position',[20,270,200,10],'FontSize',7);
hofxtext=uicontrol('Style','edit','String','0','Position',[360,90,20,15],'Visible','Off','Callback',{@ofxtext_Callback},'Units','normalized');
hofytext=uicontrol('Style','edit','String','0','Position',[380,90,20,15],'Visible','Off','Callback',{@ofytext_Callback},'Units','normalized');
betabox=uicontrol('Style','edit','String',int2str(Beta),'Position',[295,90,20,15],'Visible','On','Callback',{@setbeta_Callback},'Units','normalized');
lambdabox=uicontrol('Style','edit','String',int2str(Lambda),'Position',[315,90,20,15],'Visible','On','Callback',{@setlambda_Callback},'Units','normalized');
periodbox=uicontrol('Style','edit','String',num2str(Period,7),'Position',[335,90,20,15],'Visible','On','Callback',{@setperiod_Callback},'Units','normalized');
hpopup = uicontrol('Style','popupmenu',...
    'String',popmenu,...
    'Position',[300,50,40,25],...
    'Callback',{@popup_menu_Callback});
ha = axes('Units','Pixels','Position',[20,20,250,250]);

hOCCmenu=uicontrol('Style','popupmenu',...
    'String',OCCmenu,...
    'Position',[360,50,40,25],...
    'Callback',{@OCCmenu_Callback},'Visible','Off');
%%%%%%%%%%%%%%%Toggle selection button%%%%%%%%%%%%%%%%%%%%%%%%
 hselecttoggle=uibuttongroup('BorderType','none','Visible','off',...
                   'Position',[0.67,0.3,0.2,0.2],'SelectionChangedFcn',@bselection);
 hstv = uicontrol(hselecttoggle,'Style',...
                   'radiobutton',...
                   'String','Vertex',...
                   'Position',[0 10 70 15],...
                   'HandleVisibility','off');
%               
hstf= uicontrol(hselecttoggle,'Style','radiobutton',...
                  'String','Facet',...
                  'Position',[0 30 70 15],...
                  'HandleVisibility','on','Visible','on');
              %hselecttoggle.Visible='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the GUI.
% Change units to normalized so components resize
% automatically.
set([f,ha,hload,hloadmesh,hOCCmenu,hloadOCC,hsave,hcontour,htext,hpopup,hcheckrot,hcheckfit,hsundir,hselecttoggle,betabox,lambdabox,periodbox],...
    'Units','normalized');
cview=[1,1,1];
%Create a plot in the axes.
align([hload, hsave,hcontour,hcheckrot,hpopup],'Center','None');
[vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,1);
fig3d=trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3),shade);axis equal;xlabel('X');ylabel('Y');zlabel('Z');
hold on
axislim=axis();
TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
VN=vertexNormal(TR);
% Assign the GUI a name to appear in the window title.
set(f,'Name','Simple Mesh Editor')
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on');
set(f,'KeyPressFcn',@move_vertex)
set(f,'WindowButtonDownFcn',@VertexSelect_orig);
%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hsundir_Callback(source,eventdata)
        %Rotate view direction corresponding to sun
        axes(ha);view(rE0);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function popup_menu_Callback(source,eventdata)
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        % Set current data to the selected data set.
        cview=FT.E{val};
        [az,el]=view;
        current_view=val;
        %sprintf('Current view: %f %f\n',az,el)
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,val);
        [th1,phi1]=cart2sph(rE(1),rE(2),rE(3));
        AOtext=sprintf('Az: %3.1f El: %3.1f Phase: %3.1f Px: %3.1f\n',rad2deg(th1),rad2deg(phi1),acosd(rE*rE0'),FT.scale{val}/(1/(FT.distance{val}*149597871)*180/pi*3600));
        set(htext,'String',AOtext);
        set(fig3d,'FaceVertexCData',shade);
       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fitbutton_Callback(source,eventdata)
        %Save mesh to temporary file
        file=tempname;
        write_standard_shape_file(tlist,vlist,file);
        cdir=cd;
        cd(iniPathName);
        [status,result] = system(['adam ' iniFileName ' --checkfit ' '--shapefile ' file]);
        sprintf('%s\n',result(strfind(result,'chisq'):end))
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loadbutton_Callback(source,eventdata)
        current_view=1;
        LoadIni;
        set(hpopup,'String',popmenu,'Visible','off');
        popmenu=cell(1,nAO);
        for j=1:nAO
            popmenu{j}=strcat('AO',int2str(j));
        end
        set(hpopup,'String',popmenu,'Visible','on');
        hold off;
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,1);
        fig3d=trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3),shade);axis equal;
        hold on;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loadmeshbutton_Callback(source,eventdata)
        [FileName,PathName] = uigetfile('*.*');
        filename=strcat(PathName,FileName);
        [tlist,vlist]=read_shape(filename,1);
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,1);
        hold off
        fig3d=trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3),shade);axis equal;
        hold on
        TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
        VN=vertexNormal(TR);
        nfac=size(tlist,1);
        nvert=size(vlist,1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loadOCCbutton_Callback(source,eventdata)
        [FileName,PathName] = uigetfile('*.*','Select occultation');
        filename=strcat(PathName,FileName);
        OCC=read_OCC_struct(filename);
        nOCC=size(OCC.E,2);
        OCCOffset=zeros(1,2*nOCC);
        OCCmenu=cell(1,nOCC);
        for j=1:nOCC
            OCC.TIME{j}=OCC.TIME{j}-min_tim;
            OCCmenu{j}=strcat('OCC',int2str(j));
        end
        %  if ~ishandle(OCCFig)
        %        OCCFig=figure;
        %    end
        %    Plot_Occultation(tlist,vlist,angles,OCCoffset(2*OCCindex:2*OCCindex-1),OCC,[],OCCindex,OCCFig);
        
        OCC=modify_occ_data(OCC);
        %  hOCCmenu=uicontrol('Style','popupmenu','String',OCCmenu,...
        %            'Position',[360,50,40,25],...
        %           'Callback',{@OCCmenu_Callback},'Visible','on','Units','normalized');
        set(hOCCmenu,'String',OCCmenu,'Visible','on');
        set(hofxtext,'String','0.0','Visible','on');
        set(hofytext,'String','0.0','Visible','on');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function OCCmenu_Callback(source,eventdata)
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        OCCindex=val;
        if ~ishandle(OCCFig)
            OCCFig=figure;
        end
        set(hofxtext,'string',num2str(OCCOffset(2*OCCindex-1)));
        set(hofytext,'string',num2str(OCCOffset(2*OCCindex)));
        Plot_Occultation(tlist,vlist,angles,OCCOffset(2*OCCindex-1:2*OCCindex),OCC,[],OCCindex,OCCFig);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ofxtext_Callback(source,eventdata)
        ofx=get(source,'string');
        OCCOffset(2*OCCindex-1)=str2double(ofx);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ofytext_Callback(source,eventdata)
        ofy=get(source,'string');
        OCCOffset(2*OCCindex)=str2double(ofy);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setbeta_Callback(source,eventdata)
        Beta=str2double(get(source,'string'));
        angles(1)=deg2rad(90-str2double(get(source,'string')));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    function setlambda_Callback(source,eventdata)
        Lambda=str2double(get(source,'string'));
        angles(2)=deg2rad((str2double(get(source,'string'))));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setperiod_Callback(source,eventdata)
        Period=str2double(get(source,'string'));
        angles(3)=24*2*pi/str2double(get(source,'string'));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function savebutton_Callback(source,eventdata)
        % Display mesh plot of the currently selected data.
        [FileName,PathName] = uiputfile('*.*');
        filename=strcat(PathName,FileName);
        write_standard_shape_file(tlist,vlist,filename);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function projbutton_Callback(source,eventdata)
        %Display data and model projections
        if ~ishandle(DataFig)
            DataFig=figure;
        end
        if ~ishandle(ProjFig)
            ProjFig=figure;
            set(ProjFig,'Name','Model','NumberTitle','off');
        end
        figure(DataFig); imagesc(im{current_view});axis equal;axis xy;
        set(DataFig,'Name',Filename{current_view},'NumberTitle','off');
        ims=size(im{current_view},1);
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,current_view);
        I2=project_mesh_2d(tlist,vlist,angles,FT.TIME{current_view},FT.E{current_view},FT.E0{current_view},FT.distance{current_view},FT.up{current_view},FT.scale{current_view}/2,2*ims,[],Albedo,shade);
        figure(ProjFig); imagesc(I2);axis xy;axis equal;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%rotbutton_callback%%%%%%%%%%%%%%%%%%%
    function rotbutton_Callback(source,eventdata)
        if ~ishandle(Proj3DFig)
            Proj3DFig=figure;
            set(Proj3DFig,'Name','Model 3D projection','NumberTitle','off');
        end
        
        R=rot_matrix(angles(1),angles(2),angles(3),FT.TIME{current_view},0);
        [~,M,~,~,~]=Orthogonal_Proj(vlist,FT.E{current_view},FT.up{current_view});
        Mr=M*R';
        vlistR=(Mr*vlist')';
        clf(Proj3DFig);
        figure(Proj3DFig); Proj3DFigPlot=trisurf(tlist,vlistR(:,1),vlistR(:,2),vlistR(:,3),shade);view([0,0,1]);axis equal;axis xy;
        hold on
        %set(fig3d,'FaceVertexCData',shade);
        view([0,0,1]);
        if selected_vertex>0
            h = findobj(Proj3DFig.CurrentAxes,'Tag','pt');
            delete(h);
            figure(Proj3DFig);h = plot3(vlistR(selected_vertex,1), vlistR(selected_vertex,2), ...
                vlistR(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt');
            
            %%%%%%%%Do same for the main window
%             h = findobj(ha,'Tag','pt2');
%             delete(h);
%             set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
%                 vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
%             set(h,'Tag','pt2');
        end
        %%%%%%%%%%%%%%%%
        set(Proj3DFig, 'WindowButtonDownFcn', @VertexSelect);
        
        %set(f,'WindowButtonDownFcn',@VertexSelect_orig);
        figure(f);
%        keyboard
    end
%%%%%%%%%%%%%%%%%%%move_vertex%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function move_vertex(hObject,evt,handles)
        if ~exist('selected_vertex') || isempty(selected_vertex)
            sprintf('Select a vertex first!\n')
            return
        end
       
      % disp('key was pressed')
        if strcmp(evt.Key,'downarrow')
            
          %  disp('downarrow was pressed')
            vlist(selected_vertex,:)=vlist(selected_vertex,:)-VN(selected_vertex,:);
            TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
            VN=vertexNormal(TR);
            set(fig3d,'XData',reshape(vlist(tlist,1),nfac,3)' , 'YData',reshape(vlist(tlist,2),nfac,3)' , 'ZData',reshape(vlist(tlist,3),nfac,3)');
            refresh();
        
            refresh_projection;

            h = findobj(ha,'Tag','pt2');
            delete(h);
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
            figure(f);
        end
        if strcmp(evt.Key,'uparrow')
         %   disp('uparrow was pressed')
            
            vlist(selected_vertex,:)=vlist(selected_vertex,:)+VN(selected_vertex,:);
            TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
            VN=vertexNormal(TR);
            set(fig3d,'XData',reshape(vlist(tlist,1),nfac,3)' , 'YData',reshape(vlist(tlist,2),nfac,3)' , 'ZData',reshape(vlist(tlist,3),nfac,3)');
            refresh();
            refresh_projection;
            h = findobj(ha,'Tag','pt2');
            delete(h);
            set(0,'CurrentFigure',f);h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
           figure(f);
        end
        if strcmp(evt.Key,'leftarrow') && selected_vertex>1
            h = findobj(ha,'Tag','pt2');
            delete(h);
            selected_vertex=selected_vertex-1;
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
        end
        if strcmp(evt.Key,'rightarrow') && selected_vertex<nvert
            h = findobj(gca,'Tag','pt2');
            delete(h);
            selected_vertex=selected_vertex+1;
            h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
        end
        set(f,'WindowButtonDownFcn',@VertexSelect_orig);
       set(0,'CurrentFigure',f)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%refresh%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function refresh()
        
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,current_view);
        set(fig3d,'FaceVertexCData',shade);
    end
%%%%%%%%%%%%%%%%%%%%%%%function refresh_projection%%%%%%%%%%%%%
    function refresh_projection()
        %We call this when vlist changes
        if ~ishandle(Proj3DFig)
            return;
            end
        R=rot_matrix(angles(1),angles(2),angles(3),FT.TIME{current_view},0);
        [~,M,~,~,~]=Orthogonal_Proj(vlist,FT.E{current_view},FT.up{current_view});
        Mr=M*R';
        vlistR=(Mr*vlist')';
       
         set(Proj3DFigPlot,'XData',reshape(vlistR(tlist,1),nfac,3)' , 'YData',reshape(vlistR(tlist,2),nfac,3)' , 'ZData',reshape(vlistR(tlist,3),nfac,3)');
       
        %set(fig3d,'FaceVertexCData',shade);
       
        if selected_vertex>0
            h = findobj(Proj3DFig.CurrentAxes,'Tag','pt');
            delete(h);
           set(0,'CurrentFigure',Proj3DFig); h = plot3(vlistR(selected_vertex,1), vlistR(selected_vertex,2), ...
                vlistR(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt');
            
            %%%%%%%%Do same for the main window
%             h = findobj(ha,'Tag','pt2');
%             delete(h);
%             axes(ha); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
%                 vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
%             set(h,'Tag','pt2');
        end
        figure(f);
    end
        
%%%%%%%%%%%%%%%%%%%LoadIni%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function LoadIni()
        [iniFileName,iniPathName] = uigetfile('*.ini','Select ini');
        %keyboard
        %FileName='10_1_octdec.ini';
        cdir=cd;
        cd(iniPathName);
        %keyboard
        [tlist,vlist,E,E0,up,Date,angles,PixScale,km2arcsec,im,FT,Filename,RotAngle,min_tim,Albedo]=Parse_Ini(iniFileName,0);
         nfac=size(tlist,1);
        nvert=size(vlist,1);
        if isempty(Albedo)
            Albedo=ones(1,nfac);
        else
            fprintf('Albedo values loaded\n');
        end
        cd(cdir);
        %cd /tmp/checkfit;
       
        TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
        VN=vertexNormal(TR);
        nAO=size(FT.E,2);
        [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,1);
        current_view=1;
for j=1:nAO
    popmenu{j}=strcat('AO',int2str(j));
end
        
    Beta=90-rad2deg(angles(1));    
    Lambda=rad2deg(angles(2));
    Period=24*2*pi/angles(3);
   
    end
%%%%%%%%%%%%%%%%%%%Vertex select function%%%%%%%%%%%%%%%%%%%%%
    function VertexSelect(src, eventData)
        %Slightly modified version of
        % CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
        %   Babak Taati - May 4, 2005
        
        point = get(Proj3DFig.CurrentAxes, 'CurrentPoint'); % mouse click position
        camPos = get(Proj3DFig.CurrentAxes, 'CameraPosition'); % camera position
        camTgt = get(Proj3DFig.CurrentAxes, 'CameraTarget'); % where the camera is pointing to
        
        camDir = camPos - camTgt; % camera direction
        camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector
        
        % build an orthonormal frame based on the viewing direction and the
        % up vector (the "view frame")
        zAxis = camDir/norm(camDir);
        upAxis = camUpVect/norm(camUpVect);
        xAxis = cross(upAxis, zAxis);
        yAxis = cross(zAxis, xAxis);
        
        rot = [xAxis; yAxis; zAxis]; % view rotation
        Points=vlistR';
        
        rP = rot * Points;
        
        
        rPF = rot * point' ;
        
        % find the nearest neighbour to the clicked point
        Pindex = dsearchn(rP(1:2,:)',rPF(1:2));
        
        h = findobj(Proj3DFig.CurrentAxes,'Tag','pt'); % try to find the old point
        Pcoord= Points(:, Pindex);
        selected_vertex=Pindex;
        if isempty(h) % if it's the first click (i.e. no previous point to delete)
            
            % highlight the selected point
            set(0,'CurrentFigure',Proj3DFig); h = plot3(Pcoord(1,:), Pcoord(2,:), ...
                Pcoord(3,:), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt'); % set its Tag property for later use
            %Do the same for the main Window
            
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
            
        else % if it is not the first click
            
            delete(h); % delete the previously selected point
            
            % highlight the newly selected point
            set(0,'CurrentFigure',Proj3DFig); h = plot3(Pcoord(1,:), Pcoord(2,:), ...
                Pcoord(3,:), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt');  % set its Tag property for later use
            h = findobj(ha,'Tag','pt2');
            delete(h);
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
            
        end
       % fprintf("Selected vertex: %d\n",Pindex)
    end

%%%%%%%%%%%%%%%%%%%Vertex select function, function 2, select vertex on original figure%%%%%%%%%%%%%%%%%%%%%
    function VertexSelect_orig(src, eventData)
        %Slightly modified version of
        % CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
        %   Babak Taati - May 4, 2005
        
        point = get(ha, 'CurrentPoint'); % mouse click position
        camPos = get(ha, 'CameraPosition'); % camera position
        camTgt = get(ha, 'CameraTarget'); % where the camera is pointing to
        
        camDir = camPos - camTgt; % camera direction
        camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector
        
        % build an orthonormal frame based on the viewing direction and the
        % up vector (the "view frame")
        zAxis = camDir/norm(camDir);
        upAxis = camUpVect/norm(camUpVect);
        xAxis = cross(upAxis, zAxis);
        yAxis = cross(zAxis, xAxis);
        
        rot = [xAxis; yAxis; zAxis]; % view rotation
        Points=vlist';
        
        rP = rot * Points;
        
        
        rPF = rot * point' ;
        
        % find the nearest neighbour to the clicked point
        Pindex = dsearchn(rP(1:2,:)',rPF(1:2));
        
        h = findobj(ha,'Tag','pt2'); % try to find the old point
        Pcoord= Points(:, Pindex);
        selected_vertex=Pindex;
        if isempty(h) % if it's the first click (i.e. no previous point to delete)
            
            % highlight the selected point
             if ishandle(Proj3DFig)
            set(0,'CurrentFigure',Proj3DFig); h = plot3(vlistR(selected_vertex,1), vlistR(selected_vertex,2), ...
                vlistR(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt'); % set its Tag property for later use
             end
            %Do the same for the main Window
            
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
            
        else % if it is not the first click
            
            delete(h); % delete the previously selected point
            set(0,'CurrentFigure',f); h = plot3(vlist(selected_vertex,1), vlist(selected_vertex,2), ...
                vlist(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt2');
            % highlight the newly selected point
            if ishandle(Proj3DFig)
            h = findobj(Proj3DFig.CurrentAxes,'Tag','pt');
            delete(h);
            set(0,'CurrentFigure',Proj3DFig); h = plot3(vlistR(selected_vertex,1), vlistR(selected_vertex,2), ...
                vlistR(selected_vertex,3), 'r.', 'MarkerSize', 20);
            set(h,'Tag','pt');  % set its Tag property for later use
            end
            
            
        end
        %fprintf("Selected vertex: %d\n",Pindex)
    end
end
