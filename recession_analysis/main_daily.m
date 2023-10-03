addpath(genpath("C:\Users\flipl\dev\TOSSH\TOSSH_code"))

data = readtable("C:\Program Files (x86)\HBV-light\Data\LittleWashita_daily\data\ptq.txt");

Q = data.Q_mm_hr_; %mm/h
Q(Q<0) = 0;
t = datetime(num2str(data.Date), 'InputFormat', 'yyyyMMdd');
P = data.P_mm_hr_; 

figure('pos',[100 100 350 200])
plot(t,Q,'k-','linewidth',1.0)
xlabel('Date')
ylabel('Streamflow [mm/hr]')

figure('pos',[100 100 350 200])
plot(t,P,'k-','linewidth',1.0)
xlabel('Date')
ylabel('Rainfall [mm/hr]')

[Recession_Parameters, recession_month] =sig_RecessionAnalysis(Q, t, 'fitting_type','linear', 'plot_results', true, 'fit_individual', false); 
[Recession_Parameters, recession_month] =sig_RecessionAnalysis(Q, t, 'fitting_type','nonlinear', 'plot_results', true, 'fit_individual', true, 'recession_length', 3, 'filter_par', 2); 


[MRC_num_segments,Segment_slopes]  = sig_MRC_SlopeChanges(Q, t, 'recession_length', 6, 'filter_par',2, 'plot_results', true, 'seg_test', 2);

[Spearmans_rho]= sig_RecessionUniqueness(Q, t, 'plot_results', true) 
[Spearmans_rho]= sig_RecessionUniqueness(Q, t,'plot_results', true, 'recession_length', 3, 'filter_par', 2); 

