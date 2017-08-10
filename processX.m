clear

%% Load the pre and post calibration files
load('dpdx.mat','utau','eta','Rho')
cd('Precal'); load('summary.mat','U','V1','V2','TempK');cd ..
U_pre = U;  V1_pre= V1;  V2_pre= V2;  T_pre = TempK;
cd('Postcal'); load('summary.mat','U','V1','V2','TempK');cd ..
U_post = U; V1_post=V1;   V2_post= V2; T_post = TempK;
cd('Data');load('acquisition.mat'); cd ..
%% Compute the polynomials
poly_deg = 4;U_cutoff = 0.1;
T_ref = T_pre(1);   T_w = data.Thot;

U_all = [U_pre,U_post];
V1_all = [V1_pre,V1_post];
V2_all = [V2_pre,V2_post];
T_all = [T_pre,T_post];

cal_data = find(U_all>U_cutoff);
V_corr = @(V,Ta) V.*sqrt((T_w-T_ref)./(T_w-Ta));
f = @(P,V,T) polyval(P,V_corr(V,T));

%% V1
%Both calibrations
[P1,S1] = polyfit(V_corr(V1_all(cal_data),T_all(cal_data)),U_all(cal_data),poly_deg);
%Precalibration
[P1pre,S1pre] = polyfit(V_corr(V1_pre(U_pre>U_cutoff),T_pre(U_pre>U_cutoff)),...
    U_pre(U_pre>U_cutoff),poly_deg);
%Postcalibration
[P1post,S1post] = polyfit(V_corr(V1_post(U_post>U_cutoff),T_post(U_post>U_cutoff)),...
    U_post(U_post>U_cutoff),poly_deg);

cal_curve.P1 = P1;        cal_curve.S1 = S1;
cal_curve.P1pre = P1pre;  cal_curve.S1pre = S1pre;
cal_curve.P1post = P1post;  cal_curve.S1post = S1post;
%% V2
%Both calibrations
[P2,S2] = polyfit(V_corr(V2_all(cal_data),T_all(cal_data)),U_all(cal_data),poly_deg);
%Precalibration
[P2pre,S2pre] = polyfit(V_corr(V2_pre(U_pre>U_cutoff),T_pre(U_pre>U_cutoff)),...
    U_pre(U_pre>U_cutoff),poly_deg);
%Postcalibration
[P2post,S2post] = polyfit(V_corr(V2_post(U_post>U_cutoff),T_post(U_post>U_cutoff)),...
    U_post(U_post>U_cutoff),poly_deg);

cal_curve.P2 = P2;        cal_curve.S2 = S2;
cal_curve.P2pre = P2pre;  cal_curve.S2pre = S2pre;
cal_curve.P2post = P2post;  cal_curve.S2post = S2post;

%S.normr %=sqrt(sum((U_all(cal_data)-f(P,V_all(cal_data),T_all(cal_data))).^2))
% rsq = 1 - S.normr^2 / ((length(U_all(cal_data))-1) * var(U_all(cal_data)))...
%     .*(length(cal_data)-1)/(length(cal_data)-length(P));

%%  Plots the pre and post cals
figure(1)
clf
Vs = linspace(min(V1_all),max(V1_all),100);

%v1
plot(U_post,V_corr(V1_post,T_post),'ro')
hold on
plot(U_pre,V_corr(V1_pre,T_pre),'bo')
plot(U_all(cal_data),V_corr(V1_all(cal_data),T_all(cal_data)),'kx')
set(gca,'fontsize',24)
plot(f(P1,Vs,T_pre(1)),Vs,'k')
plot(f(P1pre,Vs,T_pre(1)),Vs,'b')
plot(f(P1post,Vs,T_pre(1)),Vs,'r')

%v2
plot(U_post,V_corr(V2_post,T_post),'rs')
plot(U_pre,V_corr(V2_pre,T_pre),'bs')
plot(U_all(cal_data),V_corr(V2_all(cal_data),T_all(cal_data)),'kv')
set(gca,'fontsize',24)
plot(f(P2,Vs,T_pre(1)),Vs,'k')
plot(f(P2pre,Vs,T_pre(1)),Vs,'b')
plot(f(P2post,Vs,T_pre(1)),Vs,'r')

xlabel('U (m/s)')
ylabel('V (Volts)')
legend('Postcal','Precal','location','southeast')
hold off
print('cal','-dpng')
%%
meanU = [data.ySet',data.ySet']*0;
varU = data.ySet'*0;
varV = data.ySet'*0;
skewU = data.ySet'*0;
var2U = data.ySet'*0;

iter = 1;
cd('Data')
for i = 1:data.numPos
    if data.y_plus(i) > 250
        fl = fopen(data.name{i},'r');
        temp = fread(fl,[data.dur*data.rate,3],'single');
        fclose(fl);
        
        F = f(P1,temp(:,2),data.TempK(i));
        G = f(P2,temp(:,3),data.TempK(i));
        
        F_fluc = F - mean(F);
        G_fluc = G - mean(G);
        iter
        uv(iter,1) = utau.^2.*(1-data.y_plus(i)./data.Re_tau);
        p(iter,1) = 1/4.*(var(F_fluc)-var(G_fluc));
        q(iter,1) = 1/4.*var(F_fluc-G_fluc);
        iter = iter + 1;
    else
        continue;
    end
end
uv
ft = fittype( 'beta*x+lambda*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

calFit = fit([uv,q],p,ft,opts);
beta =  calFit.beta;
alpha = calFit.lambda*beta;
lambda = calFit.lambda;
n_window = 2^17;

for i = 1:data.numPos
    
    fl = fopen(data.name{i},'r');
    temp = fread(fl,[data.dur*data.rate,3],'single');
    fclose(fl);
    
    F = f(P1,temp(:,2),data.TempK(i));
    G = f(P2,temp(:,3),data.TempK(i));
    meanU(i,1) = mean(F);
    meanU(i,2) = mean(G);
    F_fluc = F - mean(F);
    G_fluc = G - mean(G);
    u_fluc = (F_fluc+G_fluc)./2-(F_fluc-G_fluc).*lambda./2;
    v_fluc = (F_fluc-G_fluc)./2./beta;
    
    temp = cov(u_fluc,v_fluc);
    varU(i) = temp(1,1);
    varV(i) = temp(2,2);
    covUV(i) = temp(2,1);
    
    [phi_uu,F] = pwelch(u_fluc,n_window,[],[],data.rate);
    [phi_vv,F] = pwelch(v_fluc,n_window,[],[],data.rate);
    [phi_uv,F] = cpsd(u_fluc,v_fluc,n_window,[],[],data.rate);
    Puu(:,i) = phi_uu;
    Pvv(:,i) = phi_vv;
    Puv(:,i) = phi_uv;
end
cd ..