% a script to read and curate mahurangi data to CFE model input
% (A) Most of Mahurangi data are in mm/h --> convert to CFE input unit, m/h, by dividing the original value by 1000 
% (B) 

%% preparation
cd('G:\Shared drives\Ryoko and Hilary\SMSigxModel\analysis\1_code')
in_path = '../0_data/Mahurangi';
out_path = '../2_data_input/Mahurangi/full';

%% Read data

% project2.mat contains precip and flow data
load(fullfile(in_path, 'project2.mat'));

% Flow
flow = timetable(datetime(r6806(:,1),'ConvertFrom','datenum'), r6806(:,3)./1000); % Get the most downstream data % mm/h --> m/hr
flow.Properties.DimensionNames = {'Time'  'Variables'};
flow.Properties.VariableNames = {'Flow'};

% Precipitation
% To do: get the weight of each region
% Get the average of all the data
precip0 = mean([r6806(:,2) r6809(:,2) r6810(:,2) r6812(:,2) ...
    r6813(:,2) r6814(:,2) r6815(:,2) r6816(:,2) ...
    r6817(:,2) r6819(:,2) r6820(:,2) r6821(:,2) ...
    r6823(:,2) r6824(:,2) r6825(:,2) r6826(:,2) ...
    r6827(:,2) r6828(:,2) r6829(:,2) r6830(:,2) ...
    r6831(:,2) r6832(:,2) r6833(:,2) r6834(:,2) ...
    r6835(:,2) r6836(:,2) r6837(:,2) r6838(:,2)], 2)./1000; % mm/h --> m/h
precip = timetable(datetime(r6806(:,1),'ConvertFrom','datenum'), precip0);
precip.Properties.DimensionNames = {'time'  'Variables'};
precip.Properties.VariableNames = {'precip_rate'};

% PET
% two data almost the same. Ask Hilary but both seems to be fine.
load(fullfile(in_path, 'Mahurangi_PET.mat'));
PET = timetable(dt_nc, pet_nc(1,:)'./1000); % mm/h --> m/h
PET.Properties.DimensionNames = {'time'  'Variables'};
PET.Properties.VariableNames = {'PET'};

% Soil Moisture
load(fullfile(in_path, 'soilmoisture.mat'));

% get weight
sm_info = readtable(fullfile(in_path,'Mahurangi_site_info.xlsx'), 'Sheet', 'Soil_moisture_selection');
weight0 = sm_info.WeightBasedOnArea_calculatedByRyoko_;

weight = repelem(weight0',size(allsmup,1),1);

scale_smup = sum(weight .* ~isnan(allsmup),2);
scale_smlo = sum(weight .* ~isnan(allsmlo),2);

weighted_smup = allsmup .* weight ./ scale_smup;
weighted_smlo = allsmlo .* weight ./ scale_smlo;

% calculate weighted avg

% sm_30cm_avg_original = mean(allsmup, 2, 'omitnan');
sm_30cm_avg = sum(weighted_smup, 2, 'omitnan');
sm_30cm_avg_tt = timetable(datetime(alldates,'ConvertFrom','datenum'), sm_30cm_avg./100); % percentile --> fraction
sm_30cm_avg_tt.Properties.DimensionNames = {'Time'  'Variables'};
sm_30cm_avg_tt.Properties.VariableNames = {'Soil Moisture Content'};

% sm_60cm_avg_original = mean(allsmlo, 2, 'omitnan');
sm_60cm_avg = sum(weighted_smlo, 2, 'omitnan');
sm_60cm_avg_tt = timetable(datetime(alldates,'ConvertFrom','datenum'), sm_60cm_avg./100);  % percentile --> fraction
sm_60cm_avg_tt.Properties.DimensionNames = {'Time'  'Variables'};
sm_60cm_avg_tt.Properties.VariableNames = {'Soil Moisture Content'};

sm_avg_tt = timetable(datetime(alldates,'ConvertFrom','datenum'), mean([sm_30cm_avg sm_60cm_avg/100], 2));  % percentile --> fraction
sm_avg_tt.Properties.DimensionNames = {'Time'  'Variables'};
sm_avg_tt.Properties.VariableNames = {'VSMC'};

%% Calculate discharge data m3/s by multiplying m/h values by area and diviging it by 3600 
catchment_area = 46.65; % km2 % specific to Mahurangi
discharge0 = flow.Flow * catchment_area * 1000000.0 / 3600;
%    cfe_state.total_discharge = cfe_state.flux_Qout_m [m/h] *
%    cfe_state.catchment_area_km2 * 1000000.0 [m2] / 3600.0 [s/h]
discharge_cms = timetable(datetime(r6806(:,1),'ConvertFrom','datenum'), discharge0); % Get the most downstream data % mm/h --> m/hr
discharge_cms.Properties.DimensionNames = {'Time'  'Variables'};
discharge_cms.Properties.VariableNames = {'Total Discharge'};

%% Synchronize data
% Define the time range
flow_record = flow.Time(~isnan(flow.Flow));
precip_record = precip.time(~isnan(precip.precip_rate));
PET_record = PET.time(~isnan(PET.PET));
sm_record = sm_avg_tt.Time(~isnan(sm_avg_tt.VSMC));
% record_s = max([flow_record(1) precip_record(1) PET_record(1) sm_record(1)]);
% record_e = min([flow_record(end) precip_record(end) PET_record(end) sm_record(end)]);

% Specific to this data
record_s = datetime(1998,2,20,11,00,00);
record_e = datetime(2001,09,06,16,00,00); % for the whole record
% record_e = datetime(1998,4,20,11,00,00); % for the short experimental data
newTime = record_s:hours(1):record_e;

%% Save
% save input file
precip = retime(precip, newTime, 'linear');
PET = retime(PET, newTime, 'linear');
inputTT = synchronize(precip, PET);
inputTT.Properties.VariableNames = {'precip_rate' 'PET'};
writetimetable(inputTT, fullfile(out_path, 'mahurangi_1998_2001.csv'));

% save validation file
flow = retime(flow, newTime, 'linear');
PET = retime(PET, newTime, 'linear');
sm_avg_tt = retime(sm_avg_tt, newTime, 'linear');
discharge_cms = retime(discharge_cms, newTime, 'linear');
sm_30cm_avg_tt = retime(sm_30cm_avg_tt, newTime, 'linear');
sm_60cm_avg_tt = retime(sm_60cm_avg_tt, newTime, 'linear');
Timestep = timetable(newTime', [0:1:length(newTime)-1]');
zero_tt = timetable(newTime', zeros(length(newTime),1));

flow_cms = timetable(newTime',flow.Flow*46.65*1000000/3600); % m/h --> m3/s

valTT = synchronize(Timestep, precip, flow, sm_avg_tt, ...
                    zero_tt, zero_tt, zero_tt, zero_tt, discharge_cms);
valTT.Properties.VariableNames = {'Time Step' 'Rainfall' 'Flow' 'Soil Moisture Content' ...
    'Direct Runoff' 'GIUH Runoff' 'Lateral Flow' 'Base Flow' 'Total Discharge'};
writetimetable(valTT, fullfile(out_path, 'test_sm_basinavg.csv'));

valTT30 = synchronize(Timestep, precip, flow, sm_30cm_avg_tt, ...
                    zero_tt, zero_tt, zero_tt, zero_tt, discharge_cms);
valTT30.Properties.VariableNames = {'Time Step' 'Rainfall' 'Flow' 'Soil Moisture Content' ...
    'Direct Runoff' 'GIUH Runoff' 'Lateral Flow' 'Base Flow' 'Total Discharge'};
writetimetable(valTT30, fullfile(out_path, 'test_sm_30cmavg.csv'));

valTT60 = synchronize(Timestep, precip, flow, sm_60cm_avg_tt, ...
                    zero_tt, zero_tt, zero_tt, zero_tt, discharge_cms);
valTT60.Properties.VariableNames = {'Time Step' 'Rainfall' 'Flow' 'Soil Moisture Content' ...
    'Direct Runoff' 'GIUH Runoff' 'Lateral Flow' 'Base Flow' 'Total Discharge'};
writetimetable(valTT60, fullfile(out_path, 'test_sm_60cmavg.csv'));

% All SM data
for i = 1:size(allsmup,2)
    for k = [30 60]
        if k == 30
            smtt = timetable(datetime(alldates,'ConvertFrom','datenum'), allsmup(:,i)./100);
        elseif k ==60
            smtt = timetable(datetime(alldates,'ConvertFrom','datenum'), allsmlo(:,i)./100);
        end
        smtt = retime(smtt, newTime, 'linear');
        smtt.Properties.VariableNames = {'sm'};
        writetimetable(smtt, fullfile(out_path,'sm_data',sprintf('sm_d%02dcm_s%02d.csv',k,i))); 
    end
end

%% Visualize
title('Mahurangi observed data')
ax1 = subplot(4,1,[1 2]);
yyaxis left
plot(valTT.Time, valTT.Flow, 'DisplayName', 'Flow'); hold on;
ylabel('Flow (mm/hr)')
yyaxis right
plot(inputTT.time, inputTT.precip_rate, 'DisplayName', 'Precip'); hold on;
ylabel('Precip. (mm/hr)')
set(gca, 'YDir', 'reverse') 
ax2 = subplot(4,1,3);
plot(inputTT.time, inputTT.PET, 'DisplayName', 'PET'); hold on;
ylabel('PET (mm/hr)')
ax3 = subplot(4,1,4);
plot(valTT.Time, valTT.('Soil Moisture Content'), 'DisplayName', 'Soil moisture');
ylabel('VSWC (m^3/m^3)')
linkaxes([ax1 ax2 ax3],'x')
saveas(gcf, fullfile(out_path, 'plot', 'inputdata.png'))

