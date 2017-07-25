% Gather Daq Devices
function [Vset] = runCalibration(varargin)
%%
if nargin<1;
    [pathstr,name,ext] = fileparts(mfilename('fullpath'));
    cd(pathstr);
    direc = uigetdir;
    fname = input('Folder name: ','s');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(direc,fname);
    cd(direc);cd(fname)
end

DAQSetup
%% Select Proper Transducer
transducer = Pitot02;
ch = addDigitalChannel(daqCal,'Dev4',transducer.DChannel,'OutputOnly');% Motor Controller Voltage

%Open valve to pitot transducer
outputSingleScan(daqCal,1);
daqCal.removeChannel(length(daqCal.Channels))

ochan= MotorOut;
ichan =  {Temperature,TunnelStatic,Dantec,transducer};

%Add motor out
ch = addAnalogOutputChannel(daqCal,'Dev4',MotorOut.Channel,'Voltage');% Motor Controller Voltage
ch.Name = MotorOut.Name;
ch.Range = MotorOut.Range;

%Add input channels
for i = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,'Dev4',ichan{i}.Channel,'Voltage');% Motor Controller Voltage
    ch.Name = ichan{i}.Name;
    ch.Range = ichan{i}.Range;
end

%% Default motor to 0
daqCal.outputSingleScan(0);

%% Calibration Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daqCal.Rate = 10000;    % Data acquisition frequency
sampleDuration = 30;     % Data sample time
numPoints = 20;          % Number of samples
Vmax = 8.2;             % Max voltage (0-10V)
rampSpeed = .1;         % V/sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vout = zeros(sampleDuration*daqCal.Rate,1);
Vs = [linspace(0,2.8,10),linspace(3,Vmax,numPoints-10)];%linspace(0,Vmax,numPoints);
diffVs = [Vs(1), diff(Vs)];
Vset = 0;

%% Set the pause criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pauseTimes = Vs*0+30;       %Default wait time 30 seconds
pauseTimes(Vs <= 3) = 60; %Velocities less than ~10m/s wait 5 min
pauseTimes(Vs > 3) = 30;    %Velocities larger than ~10m/s wait 20 sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iteration
figure(1)
xlabel('U (m/s)')
ylabel('Voltage')
hold on
for i = 1:numPoints
    %Ramp
    rampTime = abs(diffVs(i))/rampSpeed;
    fprintf('Starting point - %d/%d :\n\tRamping to %0.3f for %0.3f sec \n',i,numPoints,Vs(i),rampTime)
    
    %%% Simpler Ramp
    if round(daqCal.Rate*rampTime) > 0
        queueOutputData(daqCal,linspace(Vset,Vs(i),round(daqCal.Rate*rampTime))');
        daqCal.startForeground();
    end
    Vset = Vs(i);
    fprintf('\tPausing for %d seconds \n',pauseTimes(i))
    pause(pauseTimes(i));
    fprintf('\tTaking data for %d seconds \n',sampleDuration)
    
    
    %% Data Acquisition
    queueOutputData(daqCal,Vout+Vset);
    [captured_data,time] = daqCal.startForeground();
    
    %%% Save the data
    data = struct('TempK',Temperature.cal(mean(captured_data(:,1))),...
        'Static_Pa',TunnelStatic.cal(mean(captured_data(:,2))),...
        'V',mean(captured_data(:,3)),...
        'V_std',std(captured_data(:,3)),...
        'Pitot_Pa',transducer.cal(mean(captured_data(:,4))),...
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
    
    plot(data.U,data.V,'bo')
    U(i) = data.U;
    V(i) = data.V;
    V_std(i) = data.V_std;
    TempK(i) = data.TempK;
    Static_Pa(i) = data.Static_Pa;
    Pitot_Pa(i) = data.Pitot_Pa;
    
    
end
hold off
save('all.mat','calData')
save('summary.mat','U','V','V_std','TempK','Static_Pa','Pitot_Pa','ichan')
print('cal','-dpng')
%Close valve to pitot transducer
ch = addDigitalChannel(daqCal,'Dev4',transducer.DChannel,'OutputOnly');% Motor Controller Voltage
outputSingleScan(daqCal,[Vset,0]);
daqCal.removeChannel(length(daqCal.Channels))

%PICK A VARIANCE CRITERIA TO ELMINATE TRANSITION

%DOWN RAMP
% daqCal.Rate = 100;
% rampTime = Vset/rampSpeed;
% queueOutputData(daqCal,linspace(Vset,0,round(daqCal.Rate*rampTime))');
% daqCal.startForeground();
end
