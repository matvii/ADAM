function [tlist,vlist,E,E0,up,Date,angles,PixScale,AOSize,km2arcsec]=Parse_Ini(inifile)
%Parse ini file
fd=fopen(inifile);
while ~feof(fd)
    line=fgetl(fd);
    line=strtrim(line);
       if strfind(line,'MinTim')==1
        break
    end
end
MinTim=sscanf(line,'MinTim=%f');

while ~feof(fd)
    line=fgetl(fd);
       if strfind(line,'UseAO=')==1
        break
    end
end

nAO=sscanf(line,'UseAO=%d');
Date=NaN(1,nAO);
Filename=cell(1,nAO);
PixScale=0.009942*ones(1,nAO);
AOSize=150*ones(1,nAO);
tao1='AOFile';
tao2='Date';
tao3='Pixscale';
tao4='AOSize';
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
        line=strtrim(line);
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
            AOSize(j)=sscanf(line,'AOSize=%f');
        end
        if feof(fd)
            break;
        end
    end
end
fclose(fd);
fd=fopen(inifile);
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
end


fd2=fopen(Anglefile);
Angles=fscanf(fd2,'%f %f %f');
fclose(fd2);
angles=[deg2rad(90-Angles(1)) deg2rad(Angles(2)) 2*pi/Angles(3)*24];
fclose(fd);
%Read ephm information
cao=1.731446326742403e+02;
M=dlmread(EphFile);
E0=[];
E=[];
for j=1:size(M,1)
    %If date is not set in ini file, read it from the fits file
    if isnan(Date(j))
        Date(j)=read_fits_date(Filename{j})+2400000.5;
    end
    [I,J]=min(abs(M(:,1)-Date(j)));
    if I>1e-2
        error('MJD-OBS time not found in ephm.dat')
    end
    E0(j,:)=M(J,2:4);
    E(j,:)=M(J,5:7);
    Date(j)=Date(j)-norm(E(j,:))/cao-MinTim;
end
dist=[];
km2arcsec=[];
up=zeros(j,3);
%Scaling and rot information 
for j=1:length(Date)
    E0(j,:)=E0(j,:)/norm(E0(j,:));
    dist(j)=norm(E(j,:));
    E(j,:)=E(j,:)/dist(j);
    km2arcsec(j)=1/(dist(j)*149597871)*180/pi*3600;
    up(j,:)=calculate_rotated_frame(Filename{j},E(j,:));
end
[tlist,vlist]=read_shape(ShapeFile,1);
end
