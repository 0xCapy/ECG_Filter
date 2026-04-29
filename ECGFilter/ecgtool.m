classdef ecgtool
    %ECGTOOL Small ECG helpers used by main.m.
    % I split these out because otherwise main.m becomes unreadable.

    methods (Static)

        function D = loadData(dataFile, fs)
            %LOADDATA Load the selected MAT file and force column vectors.

            if ~isfile(dataFile)
                error('Data file not found: %s\nPut ECGData1.mat in the same folder as main.m.', dataFile);
            end

            S = load(dataFile);

            if ~isfield(S, 'origSig')
                error('The MAT file does not contain origSig.');
            end

            if ~isfield(S, 'noisySig')
                error('The MAT file does not contain noisySig.');
            end

            origSig = double(S.origSig(:));
            noisySig = double(S.noisySig(:));

            if numel(origSig) ~= numel(noisySig)
                error('origSig and noisySig must have the same number of samples.');
            end

            N = numel(origSig);
            t = (0:N-1)' / fs;

            D = struct();
            D.origSig = origSig;
            D.noisySig = noisySig;
            D.N = N;
            D.fs = fs;
            D.t = t;
            D.durationSec = N / fs;
        end


        function HR = estimateHeartRate(origSig, fs)
            %ESTIMATEHEARTRATE Basic R-peak based heart-rate estimate.
            %
            % ECGData1 has a PVC-looking negative dip, so I try normal positive
            % R-peaks first and only use the inverted signal if that clearly fails.
        
            if exist('findpeaks', 'file') ~= 2
                error('findpeaks is missing. Install or activate Signal Processing Toolbox.');
            end
        
            xRaw = double(origSig(:));
            x = xRaw - median(xRaw);
        
            minPeakDistance = round(0.40 * fs);   % about 150 bpm maximum
            sigRange = max(x) - min(x);
        
            if sigRange <= 0
                error('The ECG signal has zero amplitude range.');
            end
        
            % Try positive R-peaks first.
            promList = [0.30, 0.22, 0.16, 0.12] * sigRange;
            heightList = [0.45, 0.35, 0.28, 0.20] * max(x);
        
            bestLocsPos = [];
            bestPksPos  = [];
        
            for k = 1:numel(promList)
                [pksTry, locsTry] = findpeaks(x, ...
                    'MinPeakDistance', minPeakDistance, ...
                    'MinPeakProminence', promList(k), ...
                    'MinPeakHeight', heightList(k));
        
                if numel(locsTry) >= 6
                    bestLocsPos = locsTry;
                    bestPksPos  = pksTry;
                    break;
                end
            end
        
            % Fallback only; I do not want the PVC dip to dominate the count.
            bestLocsNeg = [];
            bestPksNeg  = [];
        
            if isempty(bestLocsPos)
                for k = 1:numel(promList)
                    [pksTry, locsTry] = findpeaks(-x, ...
                        'MinPeakDistance', minPeakDistance, ...
                        'MinPeakProminence', promList(k), ...
                        'MinPeakHeight', max(0, 0.15 * max(-x)));
        
                    if numel(locsTry) >= 6
                        bestLocsNeg = locsTry;
                        bestPksNeg  = pksTry;
                        break;
                    end
                end
            end
        
            % Pick the peak set that survived the checks.
            if ~isempty(bestLocsPos)
                locs = bestLocsPos;
                pksOriginal = xRaw(locs);
                polarity = 1;
                detectedPeakValues = bestPksPos;
            elseif ~isempty(bestLocsNeg)
                locs = bestLocsNeg;
                pksOriginal = xRaw(locs);
                polarity = -1;
                detectedPeakValues = bestPksNeg;
            else
                warning('No reliable R-peaks were detected.');
                locs = [];
                pksOriginal = [];
                polarity = 1;
                detectedPeakValues = [];
            end
        
            durationSec = numel(xRaw) / fs;
            bpmCount = numel(locs) / durationSec * 60;
        
            if numel(locs) >= 2
                rrIntervals = diff(locs) / fs;
                bpmRR = 60 / mean(rrIntervals);
            else
                rrIntervals = [];
                bpmRR = NaN;
            end
        
            HR = struct();
            HR.pks = pksOriginal;
            HR.locs = locs;
            HR.times = (locs - 1) / fs;   % match t = (0:N-1)/fs
            HR.rrIntervals = rrIntervals;
            HR.bpmCount = bpmCount;
            HR.bpmRR = bpmRR;
            HR.polarity = polarity;
            HR.minPeakDistanceSamples = minPeakDistance;
            HR.detectedPeakValues = detectedPeakValues;
        end


        function F = designFilters(fs, cutoffHz, orderN, notchFreqHz, notchBandwidthHz, improvedCutoffHz)
            %DESIGNFILTERS Build the required low-pass filters and the extra improved one.
            %
            % The brief caps the order at 5.  MATLAB wants Wn against Nyquist,
            % not cycles/sample, so I keep both values in the output struct.

            if exist('fir1', 'file') ~= 2 || exist('butter', 'file') ~= 2
                error('Signal Processing Toolbox functions fir1 and butter are required.');
            end

            if fs <= 0
                error('Sampling frequency fs must be positive.');
            end

            if cutoffHz <= 0 || cutoffHz >= fs/2
                error('cutoffHz must be between 0 and the Nyquist frequency fs/2.');
            end

            if orderN < 1 || orderN > 5 || round(orderN) ~= orderN
                error('Filter order must be an integer between 1 and 5.');
            end

            if nargin < 4 || isempty(notchFreqHz)
                notchFreqHz = 60;
            end

            if nargin < 5 || isempty(notchBandwidthHz)
                notchBandwidthHz = 2;
            end

            if nargin < 6 || isempty(improvedCutoffHz)
                improvedCutoffHz = 80;
            end

            if notchFreqHz <= 0 || notchFreqHz >= fs/2
                error('notchFreqHz must be between 0 and the Nyquist frequency fs/2.');
            end

            if notchBandwidthHz <= 0 || notchBandwidthHz >= fs/2
                error('notchBandwidthHz must be positive and smaller than the Nyquist frequency.');
            end

            if improvedCutoffHz <= 0 || improvedCutoffHz >= fs/2
                error('improvedCutoffHz must be between 0 and the Nyquist frequency fs/2.');
            end

            Wn = cutoffHz / (fs / 2);
            WnImproved = improvedCutoffHz / (fs / 2);

            % Required FIR route: Hamming-window LPF.  a = 1 for FIR.
            bFIR = fir1(orderN, Wn, 'low', hamming(orderN + 1), 'scale');
            aFIR = 1;

            % Required IIR route: Butterworth, mainly for the smooth passband.
            [bIIR, aIIR] = butter(orderN, Wn, 'low');

            % My Step 11 route: notch the mains line, then use a wider 5th-order
            % Butterworth LPF.  filtfilt is applied later, so this is an
            % offline morphology check rather than just another causal filter.
            [bNotch, aNotch] = ecgtool.designNotch(fs, notchFreqHz, notchBandwidthHz);
            [bImprovedLP, aImprovedLP] = butter(orderN, WnImproved, 'low');

            F = struct();
            F.fs = fs;
            F.cutoffHz = cutoffHz;
            F.cutoffCyclesPerSample = cutoffHz / fs;
            F.Wn = Wn;
            F.maxAllowedOrder = 5;

            F.fir = struct();
            F.fir.name = 'FIR low-pass';
            F.fir.method = 'Window method with Hamming window';
            F.fir.order = orderN;
            F.fir.b = bFIR(:).';
            F.fir.a = aFIR;

            F.iir = struct();
            F.iir.name = 'IIR Butterworth low-pass';
            F.iir.method = 'Butterworth';
            F.iir.order = orderN;
            F.iir.b = bIIR(:).';
            F.iir.a = aIIR(:).';

            F.improved = struct();
            F.improved.name = 'Improved notch + zero-phase IIR low-pass';
            F.improved.method = '60 Hz notch followed by zero-phase Butterworth low-pass';
            F.improved.notchFreqHz = notchFreqHz;
            F.improved.notchBandwidthHz = notchBandwidthHz;
            F.improved.lowpassCutoffHz = improvedCutoffHz;
            F.improved.lowpassCutoffCyclesPerSample = improvedCutoffHz / fs;
            F.improved.lowpassWn = WnImproved;
            F.improved.lowpassOrder = orderN;
            F.improved.notch.b = bNotch(:).';
            F.improved.notch.a = aNotch(:).';
            F.improved.lowpass.b = bImprovedLP(:).';
            F.improved.lowpass.a = aImprovedLP(:).';
        end


        function [b, a] = designNotch(fs, f0, bandwidthHz)
            %DESIGNNOTCH Hand-made second-order notch around the mains frequency.
            %
            % Zeros sit on the unit circle at +/-f0.  r pulls the poles slightly
            % inside the circle; I scale DC gain back to 1 afterwards.

            w0 = 2 * pi * f0 / fs;
            r = 1 - pi * bandwidthHz / fs;
            r = min(max(r, 0.80), 0.9995);

            b = [1, -2 * cos(w0), 1];
            a = [1, -2 * r * cos(w0), r^2];

            dcGain = sum(b) / sum(a);
            if isfinite(dcGain) && abs(dcGain) > eps
                b = b / dcGain;
            end
        end


        function Y = applyFilters(noisySig, F)
            %APPLYFILTERS Run the noisy ECG through each filter path.

            noisySig = double(noisySig(:));

            Y = struct();
            Y.fir = filter(F.fir.b, F.fir.a, noisySig);
            Y.iir = filter(F.iir.b, F.iir.a, noisySig);

            if isfield(F, 'improved')
                if exist('filtfilt', 'file') ~= 2
                    error('filtfilt is required for the zero-phase improvement. Activate Signal Processing Toolbox.');
                end

                % Extra path for the report: 60 Hz notch first, then wider zero-phase LPF.
                % This keeps the PQRST shape much better than the basic causal IIR.
                notchOut = filtfilt(F.improved.notch.b, F.improved.notch.a, noisySig);
                Y.improved = filtfilt(F.improved.lowpass.b, F.improved.lowpass.a, notchOut);
            end
        end


        function T = filterParameterTable(F)
            %FILTERPARAMETERTABLE Filter settings table, kept for report export if needed.

            filterName = {F.fir.name; F.iir.name};
            designMethod = {F.fir.method; F.iir.method};
            order = [F.fir.order; F.iir.order];
            samplingFrequency_Hz = [F.fs; F.fs];
            cutoffFrequency_Hz = [F.cutoffHz; F.cutoffHz];
            cutoff_cycles_per_sample = [F.cutoffCyclesPerSample; F.cutoffCyclesPerSample];
            matlab_normalised_cutoff_Wn = [F.Wn; F.Wn];
            maximum_allowed_order = [F.maxAllowedOrder; F.maxAllowedOrder];

            T = table(filterName, designMethod, order, samplingFrequency_Hz, ...
                cutoffFrequency_Hz, cutoff_cycles_per_sample, ...
                matlab_normalised_cutoff_Wn, maximum_allowed_order);
        end


        function T = filterCoefficientTable(F)
            %FILTERCOEFFICIENTTABLE Put b/a coefficients into a readable table.
            %
            % FIR has only b coefficients in practice; IIR needs both b and a.

            maxLen = max([numel(F.fir.b), numel(F.iir.b), numel(F.fir.a), numel(F.iir.a)]);
            coefficientIndex = (0:maxLen-1).';

            FIR_b = NaN(maxLen, 1);
            FIR_a = NaN(maxLen, 1);
            IIR_b = NaN(maxLen, 1);
            IIR_a = NaN(maxLen, 1);

            FIR_b(1:numel(F.fir.b)) = F.fir.b(:);
            FIR_a(1:numel(F.fir.a)) = F.fir.a(:);
            IIR_b(1:numel(F.iir.b)) = F.iir.b(:);
            IIR_a(1:numel(F.iir.a)) = F.iir.a(:);

            T = table(coefficientIndex, FIR_b, FIR_a, IIR_b, IIR_a);
        end


        function T = improvementParameterTable(F)
            %IMPROVEMENTPARAMETERTABLE One-row table for the extra filter design.

            if ~isfield(F, 'improved')
                T = table();
                return;
            end

            ImprovementName = {F.improved.name};
            DesignStrategy = {F.improved.method};
            NotchFrequency_Hz = F.improved.notchFreqHz;
            NotchBandwidth_Hz = F.improved.notchBandwidthHz;
            LowpassOrder = F.improved.lowpassOrder;
            LowpassCutoff_Hz = F.improved.lowpassCutoffHz;
            LowpassCutoff_CyclesPerSample = F.improved.lowpassCutoffCyclesPerSample;
            LowpassCutoff_MATLAB_NyquistNormalised = F.improved.lowpassWn;
            ApplicationMethod = {'filtfilt forward-backward zero-phase filtering'};

            T = table(ImprovementName, DesignStrategy, NotchFrequency_Hz, ...
                NotchBandwidth_Hz, LowpassOrder, LowpassCutoff_Hz, ...
                LowpassCutoff_CyclesPerSample, ...
                LowpassCutoff_MATLAB_NyquistNormalised, ApplicationMethod);
        end


        function T = performanceMetrics(origSig, noisySig, Y, fs, HR, hfBand)
            %PERFORMANCEMETRICS Quantitative comparison used by the report.
            %
            % I kept the metrics fairly small: overall error/correlation, delay
            % around R-peaks, PQRST-window error, and how much 60 Hz/HF content
            % was reduced.  Since origSig is provided, it is used as the reference.

            if nargin < 5 || isempty(HR)
                HR = struct('locs', [], 'polarity', 1);
            end

            if nargin < 6 || isempty(hfBand)
                hfBand = [117, 142];
            end

            origSig = double(origSig(:));
            noisySig = double(noisySig(:));

            if numel(origSig) ~= numel(noisySig)
                error('origSig and noisySig must have the same length.');
            end

            signalNames = {'Noisy ECG'; 'FIR low-pass'; 'IIR Butterworth low-pass'};
            filterClass = {'Reference noisy input'; 'FIR'; 'IIR'};
            applicationMode = {'Unfiltered'; 'Causal filter'; 'Causal filter'};
            signalData = {noisySig; double(Y.fir(:)); double(Y.iir(:))};

            if isfield(Y, 'improved')
                signalNames{end+1, 1} = 'Improved notch + zero-phase IIR';
                filterClass{end+1, 1} = 'Improved';
                applicationMode{end+1, 1} = 'Zero-phase offline';
                signalData{end+1, 1} = double(Y.improved(:));
            end

            nSignals = numel(signalData);

            RMSE_mV = NaN(nSignals, 1);
            Corr = NaN(nSignals, 1);
            SNR_dB = NaN(nSignals, 1);
            SNR_gain_dB = NaN(nSignals, 1);
            RMSE_gain_pct = NaN(nSignals, 1);
            R_abs_delay_ms = NaN(nSignals, 1);
            PQRST_N = NaN(nSignals, 1);
            PQRST_RMSE_mV = NaN(nSignals, 1);
            PQRST_Corr = NaN(nSignals, 1);
            Red_60Hz_dB = NaN(nSignals, 1);
            HF_red_dB = NaN(nSignals, 1);

            Mains60HzMagnitude_mV = NaN(nSignals, 1);
            HighFreqMeanMagnitude_mV = NaN(nSignals, 1);

            for k = 1:nSignals
                x = signalData{k};
                if numel(x) ~= numel(origSig)
                    error('Signal length mismatch when computing performance metrics.');
                end

                err = x - origSig;
                RMSE_mV(k) = sqrt(mean(err.^2));
                Corr(k) = ecgtool.safeCorrelation(origSig, x);
                SNR_dB(k) = 10 * log10(sum(origSig.^2) / max(sum(err.^2), eps));

                [f, P1] = ecgtool.singleSidedSpectrum(x, fs);

                humMask = f >= 57 & f <= 63;
                if any(humMask)
                    Mains60HzMagnitude_mV(k) = max(P1(humMask));
                end

                hfMask = f >= hfBand(1) & f <= hfBand(2);
                if any(hfMask)
                    HighFreqMeanMagnitude_mV(k) = mean(P1(hfMask));
                end

                [~, R_abs_delay_ms(k), ~] = ecgtool.rPeakDelayMetrics(origSig, x, fs, HR);

                [PQRST_N(k), PQRST_RMSE_mV(k), ~, PQRST_Corr(k)] = ...
                    ecgtool.pqrstWindowMetrics(origSig, x, fs, HR, 0.25, 0.45);
            end

            rmseNoisy = RMSE_mV(1);
            snrNoisy = SNR_dB(1);
            humNoisy = Mains60HzMagnitude_mV(1);
            hfNoisy = HighFreqMeanMagnitude_mV(1);

            for k = 1:nSignals
                SNR_gain_dB(k) = SNR_dB(k) - snrNoisy;
                RMSE_gain_pct(k) = 100 * (rmseNoisy - RMSE_mV(k)) / max(rmseNoisy, eps);
                Red_60Hz_dB(k) = 20 * log10(max(humNoisy, eps) / max(Mains60HzMagnitude_mV(k), eps));
                HF_red_dB(k) = 20 * log10(max(hfNoisy, eps) / max(HighFreqMeanMagnitude_mV(k), eps));
            end

            Signal = signalNames;
            Filter_Class = filterClass;
            Application = applicationMode;

            T = table(Signal, Filter_Class, Application, ...
                RMSE_mV, Corr, SNR_dB, SNR_gain_dB, RMSE_gain_pct, R_abs_delay_ms, ...
                PQRST_N, PQRST_RMSE_mV, PQRST_Corr, ...
                Red_60Hz_dB, HF_red_dB, ...
                'VariableNames', { ...
                    'Signal', ...
                    'Filter_Class', ...
                    'Application', ...
                    'RMSE_mV', ...
                    'Corr', ...
                    'SNR_dB', ...
                    'SNR_gain_dB', ...
                    'RMSE_gain_pct', ...
                    'R_abs_delay_ms', ...
                    'PQRST_N', ...
                    'PQRST_RMSE_mV', ...
                    'PQRST_Corr', ...
                    'Red_60Hz_dB', ...
                    'HF_red_dB' ...
                });
        end


        function [f, P1] = singleSidedSpectrum(x, fs)
            %SINGLESIDEDSPECTRUM One-sided FFT magnitude for the plots and metrics.

            x = double(x(:));
            N = numel(x);
            Yfft = fft(x);
            P2 = abs(Yfft / N);
            P1 = P2(1:floor(N/2) + 1);
            if numel(P1) > 2
                P1(2:end-1) = 2 * P1(2:end-1);
            end
            f = fs * (0:floor(N/2))' / N;
        end


        function r = safeCorrelation(x, y)
            %SAFECORRELATION Pearson correlation, with zero-range protection.

            x = double(x(:));
            y = double(y(:));
            x = x - mean(x);
            y = y - mean(y);
            denom = sqrt(sum(x.^2) * sum(y.^2));
            if denom <= eps || ~isfinite(denom)
                r = NaN;
            else
                r = sum(x .* y) / denom;
            end
        end


        function [meanSignedDelayMs, meanAbsDelayMs, maxAbsDelayMs] = rPeakDelayMetrics(origSig, testSig, fs, HR)
            %RPEAKDELAYMETRICS Local R-peak timing shift against origSig.

            meanSignedDelayMs = NaN;
            meanAbsDelayMs = NaN;
            maxAbsDelayMs = NaN;

            if ~isfield(HR, 'locs') || isempty(HR.locs)
                return;
            end

            locs = HR.locs(:);
            polarity = 1;
            if isfield(HR, 'polarity') && ~isempty(HR.polarity)
                polarity = HR.polarity;
            end

            searchHalfWidth = max(2, round(0.06 * fs));
            delays = NaN(numel(locs), 1);
            N = numel(testSig);

            for i = 1:numel(locs)
                loc0 = locs(i);
                i1 = max(1, loc0 - searchHalfWidth);
                i2 = min(N, loc0 + searchHalfWidth);

                localSig = testSig(i1:i2);
                if polarity >= 0
                    [~, relIdx] = max(localSig);
                else
                    [~, relIdx] = min(localSig);
                end

                locTest = i1 + relIdx - 1;
                delays(i) = 1000 * (locTest - loc0) / fs;
            end

            delays = delays(isfinite(delays));
            if isempty(delays)
                return;
            end

            meanSignedDelayMs = mean(delays);
            meanAbsDelayMs = mean(abs(delays));
            maxAbsDelayMs = max(abs(delays));
        end



        function [nWindows, meanRMSE, meanMAE, meanCorr] = pqrstWindowMetrics(origSig, testSig, fs, HR, preSec, postSec)
            %PQRSTWINDOWMETRICS Compare shape in R-centred heartbeat windows.
            %
            % Window choice is deliberately a bit wider after the R peak so the
            % T wave is included for this dataset.

            nWindows = 0;
            meanRMSE = NaN;
            meanMAE = NaN;
            meanCorr = NaN;

            if nargin < 5 || isempty(preSec)
                preSec = 0.25;
            end
            if nargin < 6 || isempty(postSec)
                postSec = 0.45;
            end
            if ~isfield(HR, 'locs') || isempty(HR.locs)
                return;
            end

            origSig = double(origSig(:));
            testSig = double(testSig(:));
            N = numel(origSig);
            if numel(testSig) ~= N
                error('Signal length mismatch when computing PQRST-window metrics.');
            end

            preSamp = round(preSec * fs);
            postSamp = round(postSec * fs);
            locs = HR.locs(:);

            rmseVals = NaN(numel(locs), 1);
            maeVals = NaN(numel(locs), 1);
            corrVals = NaN(numel(locs), 1);

            for i = 1:numel(locs)
                centre = locs(i);
                i1 = centre - preSamp;
                i2 = centre + postSamp;

                if i1 < 1 || i2 > N || i1 >= i2
                    continue;
                end

                refWin = origSig(i1:i2);
                testWin = testSig(i1:i2);
                errWin = testWin - refWin;

                rmseVals(i) = sqrt(mean(errWin.^2));
                maeVals(i) = mean(abs(errWin));
                corrVals(i) = ecgtool.safeCorrelation(refWin, testWin);
            end

            valid = isfinite(rmseVals) & isfinite(maeVals) & isfinite(corrVals);
            nWindows = sum(valid);
            if nWindows == 0
                return;
            end

            meanRMSE = mean(rmseVals(valid));
            meanMAE = mean(maeVals(valid));
            meanCorr = mean(corrVals(valid));
        end
    end
end