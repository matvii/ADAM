function [tlist,vlist,E,E0,up,Date,angles,PixScale,km2arcsec,im,FT,RD,Filename,RotAngle,MinTim,Albedo]=Parse_Ini_with_RD(inifile,reduce)
%Parse ini file
fd=fopen(inifile);
Hapke=[];
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'MinTim')==1
        break
    end
end
MinTim=sscanf(line,'MinTim=%f');
fclose(fd);
fd=fopen(inifile);
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
    if strfind(line,'HapkeParams=')==1
        break;
    end
end
Hapke=sscanf(line,'HapkeParams=%f,%f,%f,%f,%f');
fclose(fd);
fd=fopen(inifile);
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'UseAO=')==1
        break
    end
end

nAO=sscanf(line,'UseAO=%d');
Date=NaN(1,nAO);
Filename=cell(1,nAO);
PixScale=0.009942*ones(1,nAO);
AOSizex=zeros(1,nAO);
AOSizey=zeros(1,nAO);
AORectx=zeros(1,nAO);
AORecty=zeros(1,nAO);
RotAngle=zeros(1,nAO);
FixCam=zeros(1,nAO);
im=cell(1,nAO);
tao1='AOFile';
tao2='Date';
tao3='PixScale';
tao4='AOSize';
tao5='RotAngle';
tao6='AORect';
tao=strcat('[AO','1',']');
while ~feof(fd)
    line=fgetl(fd);
       if strfind(line,tao)==1
        break
    end
    end
for j=1:nAO
    
    tao=strcat('[AO',int2str(j),']');
    taon=strcat('[AO',int2str(j+1),']');
   
    
    
    while  isempty(strfind(line,taon)==1)
    line=fgetl(fd);
    
    line=strtrim(line);
    if ~isempty(line) && line(1)=='#'
        line=fgetl(fd);
        if ~feof(fd)
        line=strtrim(line);
        else
            break;
        end
    end
       if strfind(line,tao1)==1
            Filename{j}=sscanf(line,'AOFile=%s');
       end
        if strfind(line,tao2)==1
            Date(j)=sscanf(line,'Date=%f');
        end
        if strfind(line,tao3)==1
            PixScale(j)=sscanf(line,'PixScale=%f');
        end
        if strfind(line,tao4)==1
            [temp,s]=sscanf(line,'AOSize=%f,%f');
            if s==2
            AOSizex(j)=temp(1);
            AOSizey(j)=temp(2);
            end
        end
         if strfind(line,tao6)==1
            [temp,s]=sscanf(line,'AORect=%f,%f');
            if s==2
            AORectx(j)=temp(1);
            AORecty(j)=temp(2);
            end
        end
        if strfind(line,tao5)==1
            RotAngle(j)=sscanf(line,'RotAngle=%f');
        end
        if strfind(line,'SetCamera')==1
            FixCam(j)=1;
        end
        if feof(fd)
            break;
        end
    end
end
im=cell(1,nAO);
fclose(fd);
fd=fopen(inifile);
AlbedoFile='';
Albedo='';
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'EphFile')==1
        EphFile=sscanf(line,'EphFile=%s'); 
       end
    if strfind(line,'ShapeFile')==1
        ShapeFile=sscanf(line,'ShapeFile=%s');
    end
    if strfind(line,'AnglesFile')==1
        Anglefile=sscanf(line,'AnglesFile=%s');
    end
    if strfind(line,'AlbedoFile')==1
        AlbedoFile=sscanf(line,'AlbedoFile=%s');
    end
end
fclose(fd);

fd2=fopen(Anglefile);
%keyboard
Angles=fscanf(fd2,'%f %f %f');
fclose(fd2);
angles=[deg2rad(90-Angles(1)) deg2rad(Angles(2)) 2*pi/Angles(3)*24];
%fclose(fd);
%keyboard
if ~isempty(AlbedoFile)
    try
    Albedo=dlmread(AlbedoFile);
    catch
        warning('Albedo file not found');
    end
end
%Read ephm information
cao=1.731446326742403e+02;
M=dlmread(EphFile);
E0=[];
E=[];
for j=1:nAO
    %If date is not set in ini file, read it from the fits file
    if isnan(Date(j))
        Date(j)=read_fits_date(Filename{j})+2400000.5;
       % sprintf('Date: %7.5f\n',Date(j))
    end
    [I,J]=min(abs(M(:,1)-Date(j)));
    if I>1e-2
        error('MJD-OBS time not found in ephm.dat')
    end
    E0(j,:)=M(J,2:4);
    E(j,:)=M(J,5:7);
   % sprintf('Date: %7.5f\n',Date(j)-norm(E(j,:))/cao)
    Date(j)=Date(j)-norm(E(j,:))/cao-MinTim;
    
end
dist=[];
km2arcsec=[];
up=zeros(j,3);
%Scaling and rot information 
for j=1:nAO
    E0(j,:)=E0(j,:)/norm(E0(j,:));
    dist(j)=norm(E(j,:));
    E(j,:)=E(j,:)/dist(j);
    km2arcsec(j)=1/(dist(j)*149597871)*180/pi*3600;
    if FixCam(j)==0
    up(j,:)=calculate_rotated_frame(Filename{j},E(j,:));
    else
        up(j,:)=[0,0.397748474527011,0.917494496447491];
    end
    %keyboard
    if RotAngle(j)~=0
    up(j,:)=calc_cam_angle(E(j,:),up(j,:),RotAngle(j)*pi/180);
    end
    im{j}=fitsread(Filename{j});
    
    if AORectx(j)>0 && AORecty(j)>0
        
        im{j}=im{j}(AORecty(j):AOSizey(j)-1+AORecty(j),AORectx(j):AOSizex(j)-1+AORectx(j));
    end
    if reduce==1
        im{j}=im{j}(36:115,36:115);
    end
    FT.E{j}=E(j,:);
    FT.E0{j}=E0(j,:);
    FT.TIME{j}=Date(j);
    FT.scale{j}=PixScale(j);
    FT.distance{j}=dist(j);
    FT.up{j}=up(j,:);
end
FT.HapkeParams=Hapke';
[tlist,vlist]=read_shape(ShapeFile,1);
nRD=0;
fd=fopen(inifile);
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'UseRD=')==1
        break
    end
end

nRD=sscanf(line,'UseRD=%d');
fclose(fd);
fd=fopen(inifile);
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'RDexp=')==1
        break
    end
end
RDexp=sscanf(line,'RDexp=%f');
fclose(fd);
fd=fopen(inifile);
RDDate=NaN(1,nRD);
RDFilename=cell(1,nRD);
RDPixScalef=zeros(1,nRD);
RDPixScalet=zeros(1,nRD);
RDFCx=zeros(1,nRD);
RDFCy=zeros(1,nRD);
RDSizex=zeros(1,nRD);
RDSizey=zeros(1,nRD);
RDRectx=zeros(1,nRD);
RDRecty=zeros(1,nRD);
RDFreq=zeros(1,nRD);

tao1='RDFile';
tao2='Date';
tao3='RadarFreq';
tao4='RDFromCenter';
tao=strcat('[RD','1',']');
while ~feof(fd)
    line=fgetl(fd);
       if strfind(line,tao)==1
        break
    end
    end
for j=1:nRD
    
    tao=strcat('[RD',int2str(j),']');
    taon=strcat('[RD',int2str(j+1),']');
   
    
    
    while  isempty(strfind(line,taon)==1)
    line=fgetl(fd);
    
    line=strtrim(line);
    if ~isempty(line) && line(1)=='#'
        line=fgetl(fd);
        if ~feof(fd)
        line=strtrim(line);
        else
            break;
        end
    end
       if strfind(line,tao1)==1
            RDFilename{j}=sscanf(line,'RDFile=%s');
       end
        if strfind(line,tao2)==1
            RDDate(j)=sscanf(line,'Date=%f');
        end
        if strfind(line,tao3)==1
            RDFreq(j)=sscanf(line,'RadarFreq=%f');
        end
        if strfind(line,tao4)==1
            [temp,s]=sscanf(line,'RDFromCenter=%f,%f');
            if s==2
            RDFCx(j)=temp(1);
            RDFCy(j)=temp(2);
            end
        end
         
        if feof(fd)
            break;
        end
    end
end
fclose(fd);
RDim=cell(1,nRD);
%Read ephm information
cao=1.731446326742403e+02;
M=dlmread(EphFile);

E=[];
for j=1:nRD
    %If date is not set in ini file, read it from the fits file
    if isnan(RDDate(j))
        RDDate(j)=read_fits_date_RD(RDFilename{j});
       % sprintf('Date: %7.5f\n',Date(j))
    end
    [I,J]=min(abs(M(:,1)-RDDate(j)));
    if I>1e-2
        error('MJD-OBS time not found in E')
    end
   
    E(j,:)=M(J,5:7);
   % sprintf('Date: %7.5f\n',Date(j)-norm(E(j,:))/cao)
    RDDate(j)=RDDate(j)-norm(E(j,:))/cao-MinTim;
    E(j,:)=E(j,:)/norm(E(j,:));
    
end
for j=1:nRD
    [im,~,cdelt1,cdelt2,crpix1,crpix2]=process_RD_files(RDFilename{j},RDFCx(j),RDFCy(j));
    RD.im{j}=im;
    RD.TIME{j}=RDDate(j);
    RD.E{j}=E(j,:);
    RD.df{j}=cdelt1;
    RD.dt{j}=cdelt2;
    RD.RDFreq{j}=RDFreq(j);
end
end
