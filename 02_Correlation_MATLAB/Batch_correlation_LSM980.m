% Author: Jan-Hagen Krohn, MPI for Biochemistry, 2024

%% Inputs
% List of file identifiers for Zeiss .raw files
% The specifiers R, P, and K must be stated here. NOT CHANNEL!


fileSpecifiers = [];
fileFolders = {};
namePatterns = {};
i_file = 1;
% fileFolder is folder in which all files are to be found (and where correlation functions are written as .csv files)
% namePattern is file name prefix before the R, P, K, Channel indices (without underscore!)


% This structure is used to define the input for one batch of z-scans
% The example is 4 z-scans of 12 positions each within the same GUV.
for r = 1:1
    for p = 1:48
        for k = 1:1
            fileSpecifiers = [fileSpecifiers; [r, p, k]];
            fileFolders{i_file} = '..\Data';
            namePatterns{i_file} = '20240501_104521zScanFCS_22a6cfdf448dc03f61c990bfe39fa654';
            i_file = i_file + 1;
        end
    end
end


% Channels to correlate:
% Each rows specifies a correlation operation
% Two elements per row specify the channels to correlate. 
% For example, to create autocorrelation of ch1 and
% cross-correlation between ch1 and ch2, but not ch2 autocorrelation, write:
% corrChannels = [
%   1, 1
%   1, 2
%   ];

% Note that at LSM980, channel naming is a bit more complex. Here's the
% lookup table for which index in corrChannels translates to which channel:
% 1 -> ChS1
% 2 -> ChS2
% 3 -> Ch2
% 4 -> GaAsP1 (NIR1)

corrChannels = [
    2, 2
    ];

% Is cross correlation symmetric? If yes,average forward CCF and backward
% CCF (e.g. ch1->ch2 and ch2->ch1). Note that if you choose true, it does
% not matter whether you list forward, backward or both cross-correlations
% in corrChannels. As long as at least one of them is mentioned, the
% software sorts things out without performing redundant calculations. 
crossCorrSymm = true;
% Other correlation settings
lagmin_s = 1E-6; % Minimum lag time in s
lagmax_s = 1; % Maximum lag time in s
Sampling = 12; % How many data points per factor 2 lag time span
Offset_s = 0; % Channel temporal offset in s, usually 0

nSegments = 6; % How many segment ACCs/CCCs for standard deviation calculation
Subtract_afterpulsing = true; % Subtract calibrated AP pattern, or ignore

correct_bleaching = false; % Whether to run bleaching/drift correction on photons

%% Global initialization
channels2use = unique(corrChannels);

if Subtract_afterpulsing == true
    load 'detectorsD118_980.mat'
    G_afterpulse = @(Lags,AP_char, Cntrate_Hz) ...
        (AP_char(1).*exp(-Lags./AP_char(2))+AP_char(3).*exp(-Lags./AP_char(4)))...
        ./(1+AP_char(1).*AP_char(2)+AP_char(3)*AP_char(4))./Cntrate_Hz;
end

if crossCorrSymm
    corrChannels = unique(sort(corrChannels, 2), 'rows');
end

if nSegments > 1
    disp('Standard deviation will be calculated according to Wohland method...')
else
    disp('Data will be exported without empirical standard deviation...')
end

%% File-wise processing
parfor i_file = 1:size(fileSpecifiers, 1)
    % Read-in. Note that the naming pattern of readname must be adapted to
    % the files of interest.
    curr_file_spec = fileSpecifiers(i_file,:);
    curr_name_pattern = namePatterns{i_file};
    curr_folder = fileFolders{i_file};
    readname = fullfile(curr_folder, [curr_name_pattern, '_R', num2str(curr_file_spec(1), '%0g'), '_P', num2str(curr_file_spec(2), '%0g'), '_K', num2str(curr_file_spec(3), '%0g')]);

    % Define output name based on some settings
    if correct_bleaching
        savename = [readname, '_bl'];
    else
        savename = readname;
    end
    if Subtract_afterpulsing
        savename = [savename, '_ap'];
    end
    savename = [savename, '_corr'];
    
    disp(['Processing ' [curr_name_pattern, '_R', num2str(curr_file_spec(1), '%0g'), '_P', num2str(curr_file_spec(2), '%0g'), '_K', num2str(curr_file_spec(3), '%0g')] '...'])

    % Import all channels for the chosen file.
    arrivalTimes = {};

    if correct_bleaching
        photon_weights = {};
    end % if correct_bleaching

    try
        for iChannel = 1:length(channels2use)

            switch channels2use(iChannel)
                case 1
                    channelNameSuffix = '_ChS1';
                case 2
                    channelNameSuffix = '_ChS2';
                case 3
                    channelNameSuffix = '_Ch2';
                case 4
                    channelNameSuffix = '_GaAsP1';
                otherwise
                    error('Undefined channel index. Allowed: 1->ChS1, 2->ChS2, 3->Ch2, 4->GaAsP1')
            end % switch channels2use(iChannel)

            data = readConfoCor3([readname, channelNameSuffix, '.raw']);
            arrivalTimes{iChannel} = data.ph_sync;

            if correct_bleaching
                photon_weights{iChannel} = get_blcorr_weights(data.ph_sync);
            end % if correct_bleaching

        end % for iChannel = 1:length(channels2use)

        lagmin_sync = lagmin_s .* data.TTResult_SyncRate;
        lagmax_sync = lagmax_s .* data.TTResult_SyncRate;
        Offset_sync = Offset_s .* data.TTResult_SyncRate;

        % Main loop: Perform all specified correlation operations, and export.
        for iCorr = 1:size(corrChannels, 1)
            ChNum1 = find(channels2use == corrChannels(iCorr, 1));
            ChNum2 = find(channels2use == corrChannels(iCorr, 2));

            try % Try-catch, as correlation algorithm can crash for awkward combinations of dataset and correlation settigs

                % Correlation curve (CC) calculation
                if correct_bleaching == false
                    [CC_raw, lags] = cross_corr(...
                        arrivalTimes{ChNum1}, arrivalTimes{ChNum2},...
                        lagmin_sync, lagmax_sync, Sampling, Offset_sync);

                else % correct_bleaching
                    [CC_raw, lags] = cross_corr_weights(...
                        arrivalTimes{ChNum1}, photon_weights{ChNum1}, ...
                        arrivalTimes{ChNum2}, photon_weights{ChNum2},...
                        lagmin_sync, lagmax_sync, Sampling, Offset_sync);

                end % if correct_bleaching == false

                lags = lags' ./ data.TTResult_SyncRate;

                % Offset subtraction, and afterpulsing correction if needed
                if ... % Only treat afterpulsing if...
                        Subtract_afterpulsing == true ... % ...doing so was specified by user...
                        && ChNum1 == ChNum2 .... % ...the just-calculated CC is an auto-CC...
                        && (detectors(ChNum1, 1) ~= 0 || detectors(ChNum1, 3) ~= 0) % ...and there actually is something to subtract!
                    CC = CC_raw' - 1 - ...
                        G_afterpulse(lags,detectors(ChNum1,:), numel(arrivalTimes{ChNum1})...
                        ./data.MeasDesc_AcquisitionTime.*1E3);
                elseif crossCorrSymm && ChNum1 ~= ChNum2
                    % Get backward cross-correlation and average forward and backward
                    if correct_bleaching == false
                        [CC_raw2, ~] = cross_corr(arrivalTimes{ChNum2}, arrivalTimes{ChNum1}, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                    else % correct_bleaching
                        [CC_raw2, ~] = cross_corr_weights(arrivalTimes{ChNum2}, photon_weights{ChNum2}, arrivalTimes{ChNum1}, photon_weights{ChNum1}, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                    end % if correct_bleaching == false

                    CC = (CC_raw' + CC_raw2')./2 - 1;
                else
                    % Nothing special going on: Just subtract the 1 offset.
                    CC = CC_raw' - 1;
                end % if subtractAfterpulsing == ...

                if nSegments > 1

                    % Calculation of CC standard deviation
                    segment_length = data.MeasDesc_AcquisitionTime ./ 1E3 ./ nSegments .* data.TTResult_SyncRate;
                    segmentCCs = zeros(length(CC_raw), nSegments);
                    segmentCC_m = zeros(size(CC_raw));
                    % Minimization target when amplitude-matching the segment CCs to
                    % the full-length curve 
                    fun_inner = @(x, segmentCC_m, CC_raw) sum((x .* segmentCC_m - CC_raw).^2);
                    fun_outer = @(x) fun_inner(x, segmentCC_m, CC_raw);

                    % Calculate (and scale) segment CCs
                    for m = 1:nSegments
                        Ch1_segment = arrivalTimes{ChNum1}...
                            (arrivalTimes{ChNum1} >= segment_length .* (m-1) &...
                            arrivalTimes{ChNum1} < segment_length .* (m)) - segment_length .* (m-1);
                        Ch2_segment = arrivalTimes{ChNum2}...
                            (arrivalTimes{ChNum2} >= segment_length .* (m-1) &...
                            arrivalTimes{ChNum2} < segment_length .* (m)) - segment_length .* (m-1);  

                        if correct_bleaching == false
                            [segmentCC_m, ~ ] = cross_corr(Ch1_segment,Ch2_segment, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                        else % correct_bleaching == true
                            weights_1_segment = photon_weights{ChNum1}...
                                (arrivalTimes{ChNum1} >= segment_length .* (m-1) &...
                                arrivalTimes{ChNum1} < segment_length .* (m)) - segment_length .* (m-1);
                            weights_2_segment = photon_weights{ChNum2}...
                                (arrivalTimes{ChNum2} >= segment_length .* (m-1) &...
                                arrivalTimes{ChNum2} < segment_length .* (m)) - segment_length .* (m-1);  
                            [segmentCC_m, ~ ] = cross_corr_weights(Ch1_segment, weights_1_segment, Ch2_segment, weights_2_segment, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                        end % if correct_bleaching == false

                        if crossCorrSymm && ChNum1 ~= ChNum2
                            % Backward cross-correlation
                            if correct_bleaching == false
                                [segmentCC_m2, ~ ] = cross_corr(Ch2_segment,Ch1_segment, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                            else % correct_bleaching == true
                                [segmentCC_m2, ~ ] = cross_corr_weights(Ch2_segment, weights_2_segment, Ch1_segment, weights_1_segment, lagmin_sync, lagmax_sync, Sampling, Offset_sync);
                            end % if correct_bleaching == false

                            segmentCC_m = (segmentCC_m + segmentCC_m2)./2;
                        end % if crossCorrSymm && ChNum1 ~= ChNum2
                        
                        segmentCCs(:, m) = segmentCC_m .* fminsearch(fun_outer, 1);
                    end % for iSegment = 1:nSegments

                    SD_CC = std(segmentCCs, 0, 2) ./ sqrt(nSegments);

                    % Finish up
                    if crossCorrSymm && ChNum1 ~= ChNum2
                        % Two channels used symmetrically
                        if correct_bleaching
                            cntrate_scalar = sum(photon_weights{ChNum1}) + sum(photon_weights{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        else % not correct_bleaching
                            cntrate_scalar = numel(arrivalTimes{ChNum1})+numel(arrivalTimes{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        end % if correct_bleaching

                        Cntrate = [cntrate_scalar;...
                                    cntrate_scalar;...
                                    zeros(length(lags) - 2, 1)];

                    else
                        % Only one channel used or two channels used asymmetrically
                        if correct_bleaching
                            cntrate_scalar_1 = sum(photon_weights{ChNum1}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                            cntrate_scalar_2 = sum(photon_weights{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        else % not correct_bleaching
                            cntrate_scalar_1 = numel(arrivalTimes{ChNum1}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                            cntrate_scalar_2 = numel(arrivalTimes{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        end % if correct_bleaching

                        Cntrate = [cntrate_scalar_1;...
                                    cntrate_scalar_2;...
                                    zeros(length(lags) - 2, 1)];
                    end % if crossCorrSymm...

                    Out_ChiSurf = [lags(2:end), CC(2:end), Cntrate(1:end-1), SD_CC(2:end)];

                else % if nSegments ...
                    % Finish up
                    if crossCorrSymm && ChNum1 ~= ChNum2
                        % Two channels used symmetrically
                        if correct_bleaching
                            cntrate_scalar = sum(photon_weights{ChNum1}) + sum(photon_weights{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        else % not correct_bleaching
                            cntrate_scalar = numel(arrivalTimes{ChNum1})+numel(arrivalTimes{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        end % if correct_bleaching

                        Cntrate = [cntrate_scalar;...
                                    cntrate_scalar;...
                                    zeros(length(lags) - 2, 1)];

                    else
                        % Only one channel used or two channels used asymmetrically
                        if correct_bleaching
                            cntrate_scalar_1 = sum(photon_weights{ChNum1}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                            cntrate_scalar_2 = sum(photon_weights{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        else % not correct_bleaching
                            cntrate_scalar_1 = numel(arrivalTimes{ChNum1}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                            cntrate_scalar_2 = numel(arrivalTimes{ChNum2}) ./ data.MeasDesc_AcquisitionTime .* 1E3;
                        end % if correct_bleaching

                        Cntrate = [cntrate_scalar_1;...
                                    cntrate_scalar_2;...
                                    zeros(length(lags) - 2, 1)];
                    end % if crossCorrSymm...

                    Out_ChiSurf = [lags(2:end), CC(2:end), Cntrate(1:end-1)];
                    
                end % if nSegments > 1

                writematrix(Out_ChiSurf, [savename, '_ch' num2str(channels2use(ChNum1), '%0g'), 'ch', num2str(channels2use(ChNum2), '%0g'), '.csv'])

            catch % try-catch for correlation
                disp(['Error while correlating channels ' num2str(channels2use(ChNum1), '%0g') ' and ' num2str(channels2use(ChNum2), '%0g') ' in ' readname '. Skipping channels.'])

            end % try-catch     

        end % for iCorr = 1:size(corrChannels, 1)
        
    catch % try-catch for file read-in
        disp(['Error in finding or reading ' readname '.raw. Skipping this measurement...'])
    end  % try-catch

end % for n = 1:size(files, 1)

disp('Job done.')