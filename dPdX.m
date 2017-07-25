function [ ] = dpdx(varargin)
if nargin<1;
    %Default taps to sample
    taps = 3:1:21;
else
    %user provided tap list
    taps = varargin{1};
end
%Array of the number of 'skips' the scanivalve has to make
moves = diff([1,taps]);

D = 0.1298448;  %Pipe diameter in meters
dx = 25*D/19;   %Spacing between the taps in meters
rate =10000;    %Samples rate (S/s)
dur = 10;       %Duration of sampling per tap in seconds

%Allocate memory
P = taps*0;P_std= P;TempK = P;Static_Pa = P;

DAQSetup
%Homes the scanivalve
homeScani(daqCal,ScaniHome);

%Array of the channels to sample
ichan =  {Temperature,TunnelStatic,Scanivalve};
%Add input channels
for i = 1:length(ichan)
    ch = addAnalogInputChannel(daqCal,'Dev4',ichan{i}.Channel,'Voltage');% Motor Controller Voltage
    ch.Name = ichan{i}.Name;
    ch.Range = ichan{i}.Range;
end

%Iterate over the taps
for i = 1:length(moves)
    %Skip to next tap in the array
    skipScani(daqCal,ScaniSkip,moves(i))
    pause(2)
    
    %Collect the data
    daqCal.Rate = rate;    % Data acquisition frequency
    daqCal.DurationInSeconds = dur;     % Data sample time
    [data,time] = daqCal.startForeground();
    
    %Plots the raw signal
    figure(1)
    plot(time,data)
    xlabel('Time(s)')
    ylabel('V')
    legend({'Temp','Static','Scani'})
    
    %Process the raw data
    TempK(i) = Temperature.cal(mean(data(:,1)));
    Static_Pa(i) = TunnelStatic.cal(mean(data(:,2)));
    P(i) = Scanivalve.cal(mean(data(:,3)));
    P_std(i) = std(Scanivalve.cal(data(:,3)));
    
    %Plot the linear pressure gradient
    figure(2)
    plot(taps(1:i),P(1:i),'bo-')
    xlabel('Tap#')
    ylabel('Scanivalve Pressure (Pa)')
    drawnow
end
%home the scanivalve
homeScani(daqCal,ScaniHome);

%Calculate the mean pressure gradient
DPDX1 = mean(diff(P)./(diff(taps.*dx)))
DPDX2 = fit(taps'.*dx,P','poly1');
DPDX2 = DPDX2.p1

%Determine utau and eta for the runs
if(mean(Static_Pa)<100000)
[Rho, mu] = ZSI(mean(TempK),101325);
else
[Rho, mu] = ZSI(mean(TempK),mean(Static_Pa));
end
utau = sqrt((-DPDX2./Rho)*(D./4))
eta = mu./Rho./utau*1000;

%Shove all of the data into a struct for convenience.
dpdx.P = P; dpdx.P_std = P_std;dpdx.Static_Pa = Static_Pa;
dpdx.TempK = TempK;dpdx.taps = taps;dpdx.dx = dx;
dpdx.DPDX1 = DPDX1;dpdx.DPDX2 = DPDX2;
dpdx.Rho  = Rho;dpdx.mu = mu;dpdx.eta = eta;dpdx.utau= utau;
fprintf('Re_tau = %d\n',round(D/2/eta*1000))

save('dpdx.mat','dpdx','utau','eta')
end

%Homes the scanivalve
function homeScani(daqCal,ScaniHome)
ch = addDigitalChannel(daqCal,'Dev4',ScaniHome.DChannel,'OutputOnly');
outputSingleScan(daqCal,1);
pause(2)
outputSingleScan(daqCal,0);
pause(2)
daqCal.removeChannel(length(daqCal.Channels))
end

%Skips to the next tap on the scanivalve
function skipScani(daqCal,ScaniSkip,i)
ch = addDigitalChannel(daqCal,'Dev4',ScaniSkip.DChannel,'OutputOnly');
if i > 0
    for p = 1:i
        outputSingleScan(daqCal,1);
        pause(0.5)
        outputSingleScan(daqCal,0);
        pause(0.5)
    end
end
daqCal.removeChannel(length(daqCal.Channels))
end
