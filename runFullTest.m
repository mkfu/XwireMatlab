% Check list before runing the Superpipe:
% Is the validyne demodulator on for at least 15 minutes?
% Is everthing from the traverse access port properly connected?
% Is the traverse motor power on?
% Is the power supply for pressure transducer on?
% Is the solenoid controller on?
% Where is the NSTAP in the Superpipe?
% Is the NSTAP probe connected to the Dantec?
% Has the square wave test been conducted?
% Has the singal conditioning been conducted?
% Is the Dantec in operate mode?
% Is the filter working?
% Is the pipe room terminal selected to Superpipe?
% Is the superpipe power and tower controller on?
% Does the RPM input cable reads 0 Volt?
% Is the hand held controller toggled to PC and connected?
% Is the Start button pressed?

% Okay, now you can start data acquisition in the Superpipe. 
% 
% Before turning off the Superpipe, did you move the traverse to the wall?
% 
% Questions please be directed to princetonsuperpipe@gmail.com


%% Run Motor Setup
clc
[pathstr,name,ext] = fileparts(mfilename('fullpath'));
addpath(pathstr);
cd(pathstr);

%% Testing Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.numPos =   40;        % Number of y points
data.ymin =     60e-3;      % Closest point to the wall (mm)
data.ymax =     67;          % Furthest point to the wall (mm)
data.ySet = logspace(log10(data.ymin),log10(data.ymax),data.numPos);    %Y - Location set points
data.D =        0.1298448;
data.pitot =    3.209;
data.cline =    (data.D*1000-data.ymin-data.pitot)/2;

disp('Are the following testing parameters correct [Press Enter]?')
reply = input(sprintf('y_offset: %0.3f mm\ny_max: %0.3f mm\nPitot dist: %0.3f mm\nCenterline: %0.3f mm\n',...
    data.ymin,data.ymax,data.pitot,data.cline));


data.yActual = data.ySet*0;meanU = data.ySet*0;varU = data.ySet*0;
data.TempK = [data.yActual,0];data.Static_Pa = data.TempK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Gain =     64;          %Gain on the Dantec
data.Offset =   -0.652;    %Voltage
data.R0=        128.1;
data.Rext =     171;
data.l =     	60e-3;  %mm length of wire

data.alpha =    2e-3;
data.Thot = (data.Rext/data.R0-1)./data.alpha;
disp('Are the following testing parameters correct [Press Enter]?')
reply = input(sprintf('Gain: %i\nOffset: %0.3f V\nR_0: %0.2f ohms\nRext: %0.2f ohms\n',...
    data.Gain,data.Offset,data.R0,data.Rext));


%% Sampling Parameters
data.rate =     200000;    % Data acquisition frequency
data.dur =      75;           % Data sample time sec
data.Vset =     7.2;              % Voltage Set point
rampSpeed =     .1;          % V/sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Precal file
direc = uigetdir;
fname = 'Precal';mkdir(direc,fname);
cd(direc);cd(fname)

%% Query centerline
m = traverse();
% center = input('Centerline position (mm)? : ');
% while isempty(center)
%     center = input('\n Empty input: Centerline position (mm)? : ');
% end
% [pos,~] = m.locate();
% fprintf('Current location is %0.4f mm \n\n',pos)
% fprintf('Press enter to move %0.4f mm?\n',center-pos)
% pause

% Move to centerline for Precal
[pos,~] = m.locate();
disp('Moving to Centerline')
[pos,~] = m.move(data.cline-pos);
% if abs(center - pos) > 0.3
%     [~,~] = m.move(center-pos);
% end

%% Start Precal
disp('Starting Precal')
tic
Vtemp = runCalibration(fname);
send_text_message('703-508-3338','T-Mobile','Precal Done',num2str(toc))
%%
cd(direc);cd(fname)
load('summary.mat','U','V');
poly_deg = 4;U_cutoff = 1;
[P,S] = polyfit(V(U > U_cutoff),U(U > U_cutoff),poly_deg);
cd(direc)
DAQSetup
% Ramp tunnel speed to set voltage
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
ch.Name = MotorOut.Name;
ch.Range = MotorOut.Range;
%
%Ramp to setpoint and remove channel
disp('Ramping motor down')
queueOutputData(daqCal,linspace(Vtemp,data.Vset,daqCal.Rate/rampSpeed*abs(data.Vset-Vtemp))');
daqCal.startForeground();
daqCal.removeChannel(length(daqCal.Channels));

%Wait for the speed to stabilize
disp('Pausing for 20 seconds')
pause(20);
% dPdX
disp('Finding dp/dx')
%
cd(direc)
dPdX();
cd(direc);
load('dpdx.mat', 'utau', 'eta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data folder
fname = 'Data';mkdir(direc,fname);
cd(direc);cd(fname);
% Move to wall
[pos, pos2] = m.locate();
disp(sprintf('Current location is %0.4f mm \n',pos));
m.findWall()
m.move(0.246);
%Backlash Correct()
%
for i = 1:data.numPos
    fprintf('Starting point - %d/%d :\n\tMoving to %d um\n',i,data.numPos,round(data.ySet(i)*1000))
    if i == 1
        pos = m.locate();
%         if(data.ySet(i)-pos)>0.1
%             pos = m.move(data.ySet(i)-pos);
%         end
    elseif(i > 1)
        pos = m.move(data.ySet(i)-data.ySet(i-1));
    end
    data.yActual(i) = pos+data.ymin;
    
    %Take the Temperature & Static Pressure
    ichan =  {Temperature,TunnelStatic};
    %Add input channels
    for j = 1:length(ichan)
        ch = addAnalogInputChannel(daqCal,'Dev4',ichan{j}.Channel,'Voltage');% Motor Controller Voltage
        ch.Name = ichan{j}.Name;
        ch.Range = ichan{j}.Range;
    end
    daqCal.Rate = 10000;
    daqCal.DurationInSeconds = 10;
    fprintf('\tSampling the temperature and static pressure for %d secs\n',daqCal.DurationInSeconds)
    pre_data = daqCal.startForeground();
    daqCal.removeChannel(1:length(daqCal.Channels))
    
    data.TempK(i) = Temperature.cal(mean(pre_data(:,1)));
    data.Static_Pa(i) = TunnelStatic.cal(mean(pre_data(:,2)));
    
    %Take the hotwire data
    ichan =  {Dantec};
    %Add input channels
    for j = 1:length(ichan)
        ch = addAnalogInputChannel(daqCal,'Dev4',ichan{j}.Channel,'Voltage');% Motor Controller Voltage
        ch.Name = ichan{j}.Name;
        ch.Range = ichan{j}.Range;
    end
    daqCal.Rate = data.rate;
    daqCal.DurationInSeconds = data.dur;
    
    fprintf('\tSampling the Dantec for %d secs\n',daqCal.DurationInSeconds)
    [data_hw,time] = daqCal.startForeground();
    daqCal.removeChannel(1:length(daqCal.Channels));
    
    fprintf('\tConverting Data with Precal\n')
    hwdata = Dantec.cal(P,data_hw);
    meanU(i) = mean(hwdata);
    varU(i) = var(hwdata);
    
    %Plots the raw signal
    figure(1)
    semilogx((data.yActual(1:i))./(eta),meanU(1:i)/utau,'bo-')
    xlabel('y^+')
    ylabel('U^+')
    
    %Plots the raw signal
    figure(2)
    semilogx((data.yActual(1:i))./(eta),varU(1:i)/utau^2,'bo-')
    xlabel('y^+')
    ylabel('u^{2+}')
    
    drawnow
    
    data.name{i} = sprintf('V%0.2f_Index%i_YLoc%0.2f.bin',data.Vset,i,data.ySet(i)*1000);
    fid = fopen(data.name{i},'wb');
    fwrite(fid,[time,data_hw],'single'); %
    fclose(fid);
    %fread(fopen(data.name{i},'r'),[2,inf],'ubit16');
end

%Take the Temperature & Static Pressure
ichan =  {Temperature,TunnelStatic};
%Add input channels
for j = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,'Dev4',ichan{j}.Channel,'Voltage');% Motor Controller Voltage
    ch.Name = ichan{j}.Name;
    ch.Range = ichan{j}.Range;
end
daqCal.Rate = 10000;
daqCal.DurationInSeconds = 10;
fprintf('\tSampling the temperature and static pressure for %d secs\n',daqCal.DurationInSeconds)
pre_data = daqCal.startForeground();
daqCal.removeChannel(1:length(daqCal.Channels))

data.TempK(i+1) = Temperature.cal(mean(pre_data(:,1)));
data.Static_Pa(i+1) = TunnelStatic.cal(mean(pre_data(:,2)));

%Get the average temperature during the run
data.TempK = data.TempK(1:end-1)+diff(data.TempK)./2;
data.Static_Pa = data.Static_Pa(1:end-1)+diff(data.Static_Pa)./2;
% Save the testing Data
send_text_message('703-508-3338','T-Mobile','Data taking is Done',num2str(toc))
figure(1)
    print('mean','-dpng')
figure(2)
    print('var','-dpng')

save('acquisition.mat','data')

% Ramp down for the Postcal
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
ch.Name = MotorOut.Name;
ch.Range = MotorOut.Range;
%
%Ramp to setpoint and remove channel
disp('Ramping motor down')

daqCal.Rate = 1000;
queueOutputData(daqCal,linspace(data.Vset,0,daqCal.Rate/rampSpeed*abs(data.Vset))');
daqCal.startForeground();
daqCal.removeChannel(length(daqCal.Channels));
%
cd(direc)
% Generate Postcal file
cd('Precal'); load('summary.mat','U','V');cd ..

fname = 'Postcal';mkdir(direc,fname);
cd(direc);cd(fname);
% Move to centerline for Postcal
disp('Moving to Centerline')
[pos,~] = m.move(data.cline-pos);
if abs(data.cline - pos) > 0.3
    [~,~] = m.move(data.cline-pos);
end
% Start Postcal
disp('Starting Postcal')
figure(1)
clf
plot(U,V,'rx')
hold on

Vtemp = runCalibration(fname);
send_text_message('703-508-3338','T-Mobile','Postcal is Done',num2str(toc))

% Ramp voltage down
DAQSetup
% Ramp down for the Postcal
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
ch.Name = MotorOut.Name;
ch.Range = MotorOut.Range;

disp('Ramping motor down')
queueOutputData(daqCal,linspace(Vtemp,0,daqCal.Rate/rampSpeed*abs(Vtemp))');
daqCal.startForeground();
%
cd(direc);
%release(daqCal);
% Process
%
process
send_text_message('703-508-3338','T-Mobile',sprintf('Re_tau = %d',round(data.D/2/eta*1000)),'Test Completed')
% 
% DAQSetup
% %% Ramp down for the Postcal
% %Add motor out
% DAQSetup
% disp('Adding motor channel')
% ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
% ch.Name = MotorOut.Name;
% ch.Range = MotorOut.Range;
% 
% 
% queueOutputData(daqCal,linspace(6,0,daqCal.Rate/rampSpeed*abs(6))');
% daqCal.startForeground();
% cd ..
% %release(daqCal);

%  DAQSetup
% %% Ramp down for the Postcal
% %Add motor out
% DAQSetup
% disp('Adding motor channel')
% ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
% ch.Name = MotorOut.Name;
% ch.Range = MotorOut.Range;
% 
% 
% queueOutputData(daqCal,linspace(0.6,0,daqCal.Rate/0.1*abs(0.6))');
% daqCal.startForeground();
% cd ..
%release(daqCal);