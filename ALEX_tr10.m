%% for trace viewing of 1green-1red alternating excitation
%% based on s_tr
%% Olivia Yang

%function ALEX_tr()
close all; clear all;
fclose('all');


flowtime=5;

%% Read data
pth=input('directory: ','s');
if isempty(pth)
    pth=pwd;
end
cd(pth);
save_file=pth;

fname=input('index # of filename [default=1]  ');
if isempty(fname)
    fname=1;
end
fname=num2str(fname);
disp(['hel' fname '.traces']);

timeunit=input('time unit [default=0.2 sec]  ');
if isempty(timeunit)
    timeunit=0.2;
end

%select the folder to which the files need to be saved
newfolder = [fname ' selected traces'];
mkdir(newfolder);

fid=fopen(['hel' fname '.traces'],'r');

%first line of binary file specifies length of trace
len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);

%number of traces
Ntraces=fread(fid,1,'int16');
disp('The number of traces is: ')
disp(Ntraces/2);

%raw is a linear array, looks like it
raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);

% selection=input('spots selected using frame 1-10-[1], frame 11-20-[2], or frame 1-20-[12]  ');
% if selection == 1
%     startFrame = 12;
% else
%     startFrame = 2;
% end

%alpha
LEAKAGE=0.12;

%convert into traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
Data(index)=raw(index);

fr1=[];fr2=[];t=1;
for n=1:10:len(end)-9
    if t
        fr1=[fr1 n+1:(n+8)];
    else
        fr2=[fr2 n+1:(n+8)];
    end
    t=~t;
end

rdonor=[];
racceptor=[];
gdonor=[];
gacceptor=[];
fret=[];

for i=1:(Ntraces/2)
%     rdonor(i,:)=Data(i*2-1,2:2:end);
%     racceptor(i,:)=Data(i*2,2:2:end);
%     gdonor(i,:)=Data(i*2-1,1:2:end);
%     gacceptor(i,:)=Data(i*2,1:2:end);
    
    rdonor(i,:)=Data(i*2-1,fr2);
    racceptor(i,:)=Data(i*2,fr2);
    gdonor(i,:)=Data(i*2-1,fr1);
    gacceptor(i,:)=Data(i*2,fr1);
    
    fret(i,:)=(gacceptor(i,:)-LEAKAGE*gdonor(i,:))./...
        (gacceptor(i,:)+gdonor(i,:)-LEAKAGE*gdonor(i,:));
end


%% View traces
hdl=figure;
i=0;
time=(0:length(fr1)-1)*2.5*timeunit;
%time=rtime(1:2:end-1);
cy5loss=[];
cy3loss=[];
dwelltimes=[];
steptimes={};sc=1;
while (Ntraces/2-i) > 0
    i = i+1 ;
    
    %trace window
    % red excitation
    figure(hdl);
    ax1=subplot(3,1,1);
    plot(time,rdonor(i,:),'g', time,racceptor(i,:)-LEAKAGE*rdonor(i,:),'r');
    title(['  Molecule ' num2str(i) ' of ' num2str(Ntraces/2)]);
    xlabel('Time(s)');
    ylabel('Intensity (a.u.)');
    temp=axis;temp(2)=time(end);
    axis tight;
    axis(temp);
    grid on;
    zoom on;
    
    %green excitation
    ax2=subplot(3,1,2);
    plot(time,gdonor(i,:),'g', time,gacceptor(i,:)-LEAKAGE*gdonor(i,:),'r');
    title(['green excitation']);
    xlabel('Time(s)');
    ylabel('Intensity (a.u.)');
    temp=axis;temp(2)=time(end);
    axis tight;
    axis(temp)
    grid on;
    zoom on;
    hold on;
    plot(flowtime,'v');hold off;
    
    %FRET
    ax3=subplot(3,1,3);
    plot(time,fret(i,:),'b');
    title(['FRET']);
    xlabel('Time(s)');
    ylabel('FRET efficiency');
    temp=axis;temp(2)=time(end);
    temp(3)=-0.1;
    temp(4)=1.1; 
    axis(temp);
    linkaxes([ax2,ax3],'x');
    grid on;
    zoom on;

    answer=input('press 0-back,1-go,2-save,3-lose signal times, 4-intermediate  ','s');
    
    if answer=='0'
        i=i-2;
    end

    if answer=='1'
        mol= input('which molecule do you choose:  ');
        i= mol-1;
    end
    
    %to save individual traces
    if answer=='2'
        output=[time' rdonor(i,:)' racceptor(i,:)' gdonor(i,:)' gacceptor(i,:)' fret(i,:)'];
        save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) '.dat'],'output','-ascii');
    end
    
    % for getting times of loss of green or red signal
    if answer=='3'
        disp('  left click for green, right click for red');
        [x, y, button] = ginput(2);
        if ismember(1,button) && ismember(3,button)
            for j=1:2
                if button(j)==1
                    cy3loss=[cy3loss x(j)];
                else
                    cy5loss=[cy5loss x(j)];
                end
            end
        else
            disp('   didnt save it');
        end
    end
    
    % dwell time of intermediate
    if answer=='4'
        disp('  two clicks around intermediate');
        [x, y, button] = ginput(2);
        dt=abs(x(1)-x(2));
        dwelltimes=[dwelltimes dt];
    end
    
    % save section of trace
    if answer=='5'
        disp('  two clicks');
        [x, y, button] = ginput(2);
        x(1) = round(x(1)/timeunit); % Starting point of the cut region
        x(2) = round(x(2)/timeunit); % End point of the cut region
        if x(1)<0 %prevent the script to crash when clicking ouside the trace
            x(1)=1;
        end
        if x(2)>length(time)
            x(2)=length(time);
        end
        output=[time(x(1):x(2))' gdonor(i,(x(1):x(2)))' gacceptor(i,(x(1):x(2)))'];
        save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) '.dat'],'output','-ascii');
    end
    
    %to save just fret traces
    if answer=='6'
        disp('  two clicks');
        [x, y, button] = ginput(2);
        x(1) = round(x(1)/timeunit); % Starting point of the cut region
        x(2) = round(x(2)/timeunit); % End point of the cut region
        if x(1)<0 %prevent the script to crash when clicking ouside the trace
            x(1)=1;
        end
        if x(2)>length(time)
            x(2)=length(time);
        end
        output=[time(x(1):x(2))' gdonor(i,x(1):x(2))' gacceptor(i,x(1):x(2))'];
        save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) '.dat'],'output','-ascii');
    end
    
    if answer=='7'
        disp('  click at each event, enter to end');
        [x, y, button] = ginput;
        steptimes(sc,:)={i,x'};
        sc=sc+1;
    end
    if answer=='8'
        disp('  click at each event, enter to end');
        [x, y, button] = ginput;
        steptimes_transient(sc,:)={i,x',y'};
        sc=sc+1;
    end    
    % skip to end and save things
    if answer=='9'
        break
    end
end

if ~isempty(cy3loss)
    output=[cy3loss' cy5loss'];
    save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) '_loss_times.dat'],'output','-ascii');
end

%%after things

tabl=cell2table(steptimes);
writetable(tabl,'steptimes.csv');

tabl=cell2table(steptimes);
writetable(tabl,'steptimes_transient.csv');

diff2=zeros(length(output),1);
for i=1:length(output)
    diff2(i)=output(i,1)-output(i,2);
end
% figure;
% h2=histogram(diff2,30);

count=0;zer=0;
for i=1:length(diff2)
    if diff2(i)>0
        count=count+1;
    elseif diff2(i)==0
        zer=zer+1;
    end
end
cy5first=count./length(diff2)
same=zer./length(diff2)
cy3first=1-cy5first-same

cy3loss=mean(output(:,1))-5
cy5loss=mean(output(:,2))-5

if ~isempty(dwelltimes)
    output=dwelltimes';
    save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) '_dwell_times.dat'],'output','-ascii');
end

close all;
fclose('all');


