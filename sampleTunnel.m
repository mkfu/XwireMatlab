clear
DAQXSetup
clc
pitot_zero = -0.21557;

%% Select Proper Transducer
transducer = Pitot02;

%Open valve to pitot transducer
ch = addDigitalChannel(daqCal,transducer.Ddev,transducer.DChannel,'OutputOnly');% Motor Controller Voltage
outputSingleScan(daqCal,1);
daqCal.removeChannel(length(daqCal.Channels))
pause(1)

%Input channels
ichan =  {Temperature,TunnelStatic,transducer,hw1,hw2};
%Add input channels
for i = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,ichan{i}.dev,ichan{i}.Channel,'Voltage');% Motor Controller Voltage
    ch.Name = ichan{i}.Name;
    ch.Range = ichan{i}.Range;
end
daqCal.Rate = 10000;
daqCal.DurationInSeconds = 3;

[captured_data,time] = daqCal.startForeground();
data = struct('TempK',Temperature.cal(mean(captured_data(:,1))),...
    'Static_Pa',TunnelStatic.cal(mean(captured_data(:,2))),...
    'Pitot_Pa',transducer.cal(mean(captured_data(:,3))),...
    'V1',mean(captured_data(:,4)),...
    'V1_std',std(captured_data(:,4)),...
    'V2',mean(captured_data(:,5)),...
    'V2_std',std(captured_data(:,5)));

if(mean(data.Static_Pa)<100000)
    [Rho, mu] = ZSI(mean(data.TempK),101325);
else
    [Rho, mu] = ZSI(mean(data.TempK),mean(data.Static_Pa));
end
data.rho = Rho;
data.Mu = mu;
data.U = sqrt(2/data.rho*(data.Pitot_Pa - transducer.cal(pitot_zero)));

daqCal.removeChannel(1:length(daqCal.Channels))


%%
fprintf('Temp (K): %0.2f \nStatic (Pa): %0.2f \nPitot (Pa): %0.2f \nU (m/s): %0.2f\nV1: %0.2f\nV2: %0.2f\n',...
    data.TempK,data.Static_Pa,data.Pitot_Pa,data.U,data.V1,data.V2);
fprintf('Rho (kg/m^3): %0.2f\nMu (Pa s): %d \n',data.rho,data.Mu)
fprintf('Pitot (V): %0.5f\n',mean(captured_data(:,3)))
fprintf('Re_D = %i\n', 0.1298448*data.U./1.2.*data.rho./data.Mu);
ch = addDigitalChannel(daqCal,transducer.Ddev,transducer.DChannel,'OutputOnly');% Motor Controller Voltage
%Close valve to pitot transducer
outputSingleScan(daqCal,0);
daqCal.removeChannel(1:length(daqCal.Channels))


