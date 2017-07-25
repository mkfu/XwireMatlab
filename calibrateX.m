% Gather Daq Devices
function [Vset] = runCalibration(varargin)
%%
rampSpeed = .1;         % V/sec
DAQXSetup
%% Select Proper Transducer
transducer = Pitot02;

%Open valve to pitot transducer
ch = addDigitalChannel(daqCal,transducer.Ddev,transducer.DChannel,'OutputOnly');% Motor Controller Voltage
outputSingleScan(daqCal,1);
daqCal.removeChannel(length(daqCal.Channels))

%Motor Out
ochan= MotorOut;
ch = addAnalogOutputChannel(daqCal,MotorOut.dev,MotorOut.Channel,'Voltage');% Motor Controller Voltage
ch.Name = MotorOut.Name;
ch.Range = MotorOut.Range;

%Input channels
ichan =  {Temperature,TunnelStatic,transducer,hw1,hw2};

%Add input channels
for i = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,ichan{i}.dev,ichan{i}.Channel,'Voltage');% Motor Controller Voltage
    ch.Name = ichan{i}.Name;
    ch.Range = ichan{i}.Range;
end

%% Default motor to 0
daqCal.outputSingleScan(0);

%% Calibration Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daqCal.Rate = 10000;    % Data acquisition frequency

if nargin==0;
    [pathstr,name,ext] = fileparts(mfilename('fullpath'));
    cd(pathstr);
    direc = uigetdir;
    fname = input('Folder name: ','s');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(direc,fname);
    cd(direc);cd(fname)
    
    %Default Calibration Parameters
    Vmax = 8.2;             % Max voltage (0-10V)
    numPoints = 20;          % Number of samples
        
    calSet.sampleDuration = 30;     % Data sample time
    calSet.Vs = [linspace(0,2.8,10),linspace(3,Vmax,numPoints-10)];%linspace(0,Vmax,numPoints);
    
elseif nargin==1;
    [pathstr,name,ext] = fileparts(mfilename('fullpath'));
    cd(pathstr);
    direc = uigetdir;
    fname = input('Folder name: ','s');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(direc,fname);
    cd(direc);cd(fname)
    calSet = varargin{1};
else
    calSet = varargin{1};
end


% sampleDuration = 30;     % Data sample time
% numPoints = 20;          % Number of samples
% Vmax = 8.2;             % Max voltage (0-10V)
% rampSpeed = .1;         % V/sec
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vout = zeros(calSet.sampleDuration*daqCal.Rate,1);
diffVs = [calSet.Vs(1), diff(calSet.Vs)];
Vset = 0;

%% Set the pause criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pauseTimes = calSet.Vs*0+30;       %Default wait time 30 seconds
% pauseTimes(calSet.Vs <= 3) = 60; %Velocities less than ~10m/s wait 5 min
% pauseTimes(calSet.Vs > 3) = 30;    %Velocities larger than ~10m/s wait 20 sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iteration
figure(1)
xlabel('U (m/s)')
ylabel('Voltage')
hold on
for i = 1:numPoints
    %Ramp
    rampTime = abs(diffVs(i))/rampSpeed;
    fprintf('Starting point - %d/%d :\n\tRamping to %0.3f for %0.3f sec \n',i,numPoints,calSet.Vs(i),rampTime)
    
    %%% Simpler Ramp
    if round(daqCal.Rate*rampTime) > 0
        queueOutputData(daqCal,linspace(Vset,calSet.Vs(i),round(daqCal.Rate*rampTime))');
        daqCal.startForeground();
    end
    Vset = calSet.Vs(i);
    fprintf('\tPausing for %d seconds \n',pauseTimes(i))
    pause(pauseTimes(i));
    fprintf('\tTaking data for %d seconds \n',calSet.sampleDuration)
    
    
    %% Data Acquisition
    queueOutputData(daqCal,Vout+Vset);
    [captured_data,time] = daqCal.startForeground();
    
    %%% Save the data
    data = struct('TempK',Temperature.cal(mean(captured_data(:,1))),...
        'Static_Pa',TunnelStatic.cal(mean(captured_data(:,2))),...
        'V1',mean(captured_data(:,4)),...
        'V1_std',std(captured_data(:,4)),...
        'V2',mean(captured_data(:,5)),...
        'V2_std',std(captured_data(:,5)),...
        'Pitot_Pa',transducer.cal(mean(captured_data(:,3))),...
        'Raw',captured_data,...
        'Rate',daqCal.Rate,...
        'sampleDuration',sampleDuration);
    if(mean(data.Static_Pa)<100000)
        [Rho, mu] = ZSI(mean(data.TempK),101325);
    else
        [Rho, mu] = ZSI(mean(data.TempK),mean(data.Static_Pa));
    end
    data.rho = Rho;
    if i == 1
        calData{1} = data;
    end
    %data.Pitot_Pa = data.Pitot_Pa;% - calData{1}.Pitot_Pa;
    data.U = sqrt(2/data.rho*(data.Pitot_Pa - calData{1}.Pitot_Pa));
    calData{i} = data;
    tempName= sprintf('Raw%d.mat',i);
    fprintf('\tSaving Data as %s \n\n',tempName)
    save(tempName,'data');
    
    plot(data.U,data.V1,'bo')
    hold on
    plot(data.U,data.V2,'ro')
    hold off
    
    U(i) = data.U;
    V1(i) = data.V1;    V1_std(i) = data.V1_std;
    V2(i) = data.V2;    V2_std(i) = data.V2_std;
    
    TempK(i) = data.TempK;
    Static_Pa(i) = data.Static_Pa;
    Pitot_Pa(i) = data.Pitot_Pa;
    
    
end
hold off
save('all.mat','calData')
save('summary.mat','U','V','V_std','TempK','Static_Pa','Pitot_Pa','ichan')
print('cal','-dpng')
%Close valve to pitot transducer
outputSingleScan(daqCal,Vset);
daqCal.removeChannel(1:length(daqCal.Channels))

ch = addDigitalChannel(daqCal,transducer.Ddev,transducer.DChannel,'OutputOnly');% Motor Controller Voltage
%Close valve to pitot transducer
outputSingleScan(daqCal,0);
daqCal.removeChannel(1:length(daqCal.Channels))

%PICK A VARIANCE CRITERIA TO ELMINATE TRANSITION

%DOWN RAMP
% daqCal.Rate = 100;
% rampTime = Vset/rampSpeed;
% queueOutputData(daqCal,linspace(Vset,0,round(daqCal.Rate*rampTime))');
% daqCal.startForeground();
end
