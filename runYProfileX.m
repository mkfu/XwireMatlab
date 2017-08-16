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
data.numPos =   20;        % Number of y points
data.ymin =     1.783;      % Closest point to the wall (mm)
data.ymax =     67;          % Furthest point to the wall (mm)
data.ySet =     logspace(log10(data.ymin),log10(data.ymax),data.numPos);    %Y - Location set points
data.D =        0.1298448;
data.pitot =    6.33;      % Pitot center distance to the wall (mm)
data.cline =    (data.D*1000-data.ymin-data.pitot)/2;

disp('Are the following testing parameters correct [Press Enter]?')
reply = input(sprintf('y_offset: %0.3f mm\ny_max: %0.3f mm\nPitot dist: %0.3f mm\nCenterline: %0.3f mm\n',...
    data.ymin,data.ymax,data.pitot,data.cline));

%Pre-allocated memory
data.yActual = data.ySet*0;meanU = data.ySet*0;varU = data.ySet*0;
data.TempK = [data.yActual,0];data.Static_Pa = data.TempK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data.Gain =     64;          %Gain on the Dantec
%data.Offset =   -0.602;    %Voltage

data.R0_1=        108.6;
data.Rext_1 =     134.2;
data.Gain_1 =     32;          %Gain on the Dantec
data.Offset_1 =   -0.734;    %Voltage

data.R0_2=        105.2;
data.Rext_2 =     134.2;
data.Gain_2 =     32;          %Gain on the Dantec
data.Offset_2 =   -0.735;    %Voltage

data.l =     	60e-3;  %mm length of wire

data.alpha =    2e-3;

data.Thot_1 = (data.Rext_1/data.R0_1-1)./data.alpha;
data.Thot_2 = (data.Rext_2/data.R0_2-1)./data.alpha;

disp('Are the following testing parameters correct [Press Enter]?')
reply = input(sprintf('Gain 1: %i\nOffset 1: %0.3f V\nGain 2: %i\nOffset 2: %0.3f V\nR_0_1: %0.2f ohms\nRext_1: %0.2f ohms\nR_0_2: %0.2f ohms\nRext_2: %0.2f ohms\n',...
    data.Gain_1,data.Offset_1,data.Gain_2,data.Offset_2,data.R0_1,data.Rext_1,data.R0_2,data.Rext_2));


%% Sampling Parameters
data.rate =     300000;    % Data acquisition frequency
data.dur =      15;        % Data sample time sec
data.Vset =     3.2;       % Voltage Set point
rampSpeed =     .1;        % V/sec
Vmax =          4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calSet.sampleDuration = 10;     % Data sample time
%calSet.Vs = [linspace(0,2,10),linspace(4,Vmax,6)];%linspace(0,Vmax,numPoints);
calSet.Vs = linspace(0,Vmax,18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Precal file
direc = uigetdir;
fname = 'Precal';mkdir(direc,fname);
cd(direc);cd(fname)

%% Query centerline
m = traverse();

% Move to centerline for Precal
[pos,~] = m.locate();
disp('Moving to Centerline')
[pos,~] = m.move(data.cline-pos);
if abs(data.cline - pos) > 0.3
    [~,~] = m.move(data.cline-pos);
end

%% Start Precal
disp('Starting Precal')
tic
Vtemp = calibrateX(calSet,fname);
%%
cd(direc);cd(fname)
% load('summary.mat','U','V');
% poly_deg = 4;U_cutoff = 1;
% [P,S] = polyfit(V(U > U_cutoff),U(U > U_cutoff),poly_deg);
cd(direc)
DAQXSetup
% Ramp tunnel speed to set voltage
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,MotorOut.dev,MotorOut.Channel,'Voltage');% Motor Controller Voltage
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
%%
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
        ch = addAnalogInputChannel(daqCal,ichan{j}.dev,ichan{j}.Channel,'Voltage');% Motor Controller Voltage
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
    ichan =  {hw1,hw2};
    %Add input channels
    for j = 1:length(ichan)
        ch = addAnalogInputChannel(daqCal,ichan{j}.dev,ichan{j}.Channel,'Voltage');% Motor Controller Voltage
        ch.Name = ichan{j}.Name;
        ch.Range = ichan{j}.Range;
    end
    daqCal.Rate = data.rate;
    daqCal.DurationInSeconds = data.dur;
    
    fprintf('\tSampling the Dantec for %d secs\n',daqCal.DurationInSeconds)
    [data_hw,time] = daqCal.startForeground();
    fprintf('\tSampling the Dantec for %d secs\n',daqCal.DurationInSeconds)
    [data_hw2,time2] = daqCal.startForeground();
    daqCal.removeChannel(1:length(daqCal.Channels));
    data.name{i} = sprintf('V%0.2f_Index%i_YLoc%0.2f.bin',data.Vset,i,data.ySet(i)*1000);
    fid = fopen(data.name{i},'wb');
    fwrite(fid,[time,data_hw],'single');
    fwrite(fid,[time2,data_hw2],'single');%
    fclose(fid);
%     fprintf('\tConverting Data with Precal\n')
%     hwdata = Dantec.cal(P,data_hw);
%     meanU(i) = mean(hwdata);
%     varU(i) = var(hwdata);
%     
%     Plots the raw signal
%     figure(1)
%     semilogx((data.yActual(1:i))./(eta),meanU(1:i)/utau,'bo-')
%     xlabel('y^+')
%     ylabel('U^+')
%     
%     Plots the raw signal
%     figure(2)
%     semilogx((data.yActual(1:i))./(eta),varU(1:i)/utau^2,'bo-')
%     xlabel('y^+')
%     ylabel('u^{2+}')
%     
%     drawnow
    

    %fread(fopen(data.name{i},'r'),[3,inf],'ubit16');
end

%Take the Temperature & Static Pressure
ichan =  {Temperature,TunnelStatic};
%Add input channels
for j = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,ichan{j}.dev,ichan{j}.Channel,'Voltage');% Motor Controller Voltage
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
data.y_plus = data.yActual./eta;
data.Re_tau = data.D/2/eta*1000;

%Get the average temperature during the run
data.TempK = data.TempK(1:end-1)+diff(data.TempK)./2;
data.Static_Pa = data.Static_Pa(1:end-1)+diff(data.Static_Pa)./2;
% Save the testing Data
% figure(1)
% print('mean','-dpng')
% figure(2)
% print('var','-dpng')

save('acquisition.mat','data')

% Ramp down for the Postcal
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,MotorOut.dev,MotorOut.Channel,'Voltage');% Motor Controller Voltage
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
% figure(1)
% clf
% plot(U,V,'rx')
% hold on

Vtemp = calibrateX(calSet,fname);

% Ramp voltage down
DAQXSetup
% Ramp down for the Postcal
%Add motor out
disp('Adding motor channel')
ch = addAnalogOutputChannel(daqCal,MotorOut.dev,MotorOut.Channel,'Voltage');% Motor Controller Voltage
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
processX

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
