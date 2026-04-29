classdef ecgplot


    methods (Static)

        function fig = timeDomain(t, origSig, noisySig, ~)
            %TIMEDOMAIN Full 10 s view used at the start of the report.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);
            noisySig = noisySig(:);

            fig = ecgplot.newFigure('Time-Domain ECG Signals', 'twoPanel');
            xLimFull = ecgplot.timeLimits(t);
            yLimCommon = ecgplot.calcYLim([origSig; noisySig], 0.13);

            ax1 = axes(fig, 'Position', st.axPos.upperPanel);
            plot(ax1, t, origSig, ...
                'Color', st.col.original, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);
            ecgplot.formatAxes(ax1);
            ecgplot.panelTitle(ax1, 'Original ECG Signal');
            ylabel(ax1, 'Voltage (mV)', st.labelArgs{:});
            xlabel(ax1, 'Time (s)', st.labelArgs{:});
            xlim(ax1, xLimFull);
            ylim(ax1, yLimCommon);
            ecgplot.setSecondTicks(ax1, xLimFull, 1.0);

            ax2 = axes(fig, 'Position', st.axPos.lowerPanel);
            plot(ax2, t, noisySig, ...
                'Color', st.col.noisy, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.noisy);
            ecgplot.formatAxes(ax2);
            ecgplot.panelTitle(ax2, 'Noisy ECG Signal');
            ylabel(ax2, 'Voltage (mV)', st.labelArgs{:});
            xlabel(ax2, 'Time (s)', st.labelArgs{:});
            xlim(ax2, xLimFull);
            ylim(ax2, yLimCommon);
            ecgplot.setSecondTicks(ax2, xLimFull, 1.0);

            linkaxes([ax1, ax2], 'x');
        end


        function fig = timeDomainZoom(t, origSig, noisySig, zoomRange, ~)
            %TIMEDOMAINZOOM Shorter view so the noisy overlay is actually readable.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);
            noisySig = noisySig(:);
            idx = t >= zoomRange(1) & t <= zoomRange(2);

            fig = ecgplot.newFigure('Time-Domain ECG Comparison', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            hNoisy = plot(ax, t(idx), noisySig(idx), ...
                'Color', st.col.noisy, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.noisy);

            hOrig = plot(ax, t(idx), origSig(idx), ...
                'Color', st.col.originalSoft, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.zoomOriginal);

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, sprintf('ECG Comparison, %.1f--%.1f s', zoomRange(1), zoomRange(2)));
            xlabel(ax, 'Time (s)', st.labelArgs{:});
            ylabel(ax, 'Voltage (mV)', st.labelArgs{:});
            xlim(ax, zoomRange);
            ylim(ax, ecgplot.calcYLim([origSig(idx); noisySig(idx)], 0.16));
            ecgplot.setSecondTicks(ax, zoomRange, 0.25);

            lgd = legend(ax, [hOrig, hNoisy], {'Original ECG', 'Noisy ECG'}, ...
                'Location', 'southeast', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);
        end


        function fig = rPeakDetection(t, origSig, HR, ~)
            %RPEAKDETECTION Show the peaks used for the heart-rate estimate.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);

            fig = ecgplot.newFigure('R-Peak Detection', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            hSig = plot(ax, t, origSig, ...
                'Color', st.col.original, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);

            if ~isempty(HR.locs)
                hPk = plot(ax, HR.times, HR.pks, 'o', ...
                    'LineStyle', 'none', ...
                    'MarkerSize', st.marker.size, ...
                    'LineWidth', st.marker.lineWidth, ...
                    'MarkerEdgeColor', st.col.peak, ...
                    'MarkerFaceColor', 'none');

                % Tiny legend trick: this adds the bpm value without drawing another line.
                hNote = plot(ax, NaN, NaN, ...
                    'LineStyle', 'none', ...
                    'Marker', 'none');
            else
                hPk = [];
                hNote = [];
            end

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, 'R-Peak Detection and Heart-Rate Estimation');
            xlabel(ax, 'Time (s)', st.labelArgs{:});
            ylabel(ax, 'Voltage (mV)', st.labelArgs{:});
            xlim(ax, ecgplot.timeLimits(t));
            ylim(ax, ecgplot.calcYLim(origSig, 0.18));
            ecgplot.setSecondTicks(ax, ecgplot.timeLimits(t), 1.0);

            if ~isempty(HR.locs)
                lgd = legend(ax, [hSig, hPk, hNote], ...
                    {'Original ECG', 'Detected R-peaks', sprintf('Heart rate = %.1f bpm', HR.bpmCount)}, ...
                    'Location', 'southeast', ...
                    'Box', 'on', ...
                    'Interpreter', 'none');
            else
                lgd = legend(ax, hSig, {'Original ECG'}, ...
                    'Location', 'southeast', ...
                    'Box', 'on', ...
                    'Interpreter', 'none');
            end
            ecgplot.formatLegend(lgd);
        end


        function [f, P1] = singleSidedSpectrum(x, fs)
            %SINGLESIDEDSPECTRUM One-sided spectrum, same scaling as in the report.

            x = double(x(:));
            N = numel(x);

            Y = fft(x);
            P2 = abs(Y / N);
            P1 = P2(1:floor(N/2) + 1);

            if numel(P1) > 2
                P1(2:end-1) = 2 * P1(2:end-1);
            end

            f = fs * (0:floor(N/2))' / N;
        end


        function fig = dftNoisySignal(noisySig, fs, ~)
            %DFTNOISYSIGNAL Basic DFT plot kept for checking.

            st = ecgplot.style();
            [f, P1] = ecgplot.singleSidedSpectrum(noisySig, fs);

            yMax = max(P1);
            if ~isfinite(yMax) || yMax <= 0
                yMax = 1;
            end

            yPadTop = 0.08 * yMax;
            yPadBottom = 0.03 * yMax;

            fig = ecgplot.newFigure('DFT of Noisy ECG Signal', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);

            plot(ax, f, P1, ...
                'Color', st.col.noisy, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, 'DFT Magnitude Spectrum of Noisy ECG Signal');
            xlabel(ax, 'Frequency (Hz)', st.labelArgs{:});
            ylabel(ax, 'Magnitude (mV)', st.labelArgs{:});

            xlim(ax, [0, fs/2]);
            ylim(ax, [-yPadBottom, yMax + yPadTop]);

            freqTicks = 0:20:(fs/2);
            if numel(freqTicks) >= 2
                set(ax, 'XTick', freqTicks);
            end
        end


        function fig = dftAnnotated(noisySig, fs, cutoffHz, ~)


            st = ecgplot.style();
            [f, P1] = ecgplot.singleSidedSpectrum(noisySig, fs);

            plotMask = f <= 150;
            fPlot = f(plotMask);
            PPlot = P1(plotMask);

            yMax = max(PPlot);
            if ~isfinite(yMax) || yMax <= 0
                yMax = 1;
            end

            yPadTop = 0.13 * yMax;
            yPadBottom = 0.03 * yMax;
            yLow = -yPadBottom;
            yHigh = yMax + yPadTop;
            yRange = yHigh - yLow;

            % Find the local mains-hum peak near 60 Hz rather than assuming exactly 60.0.
            humMask = f >= 57 & f <= 63;
            humFreq = 60;
            humMag = interp1(f, P1, humFreq, 'linear', 'extrap');
            if any(humMask)
                fHum = f(humMask);
                PHum = P1(humMask);
                [humMag, idxHum] = max(PHum);
                humFreq = fHum(idxHum);
            end

            % HF noise band I observed in this dataset.
            hfBand = [117, 142];

            fig = ecgplot.newFigure('Annotated DFT of Noisy ECG Signal', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            % Shade first so it stays behind the spectrum.
            hHF = patch(ax, ...
                [hfBand(1), hfBand(2), hfBand(2), hfBand(1)], ...
                [yLow, yLow, yHigh, yHigh], ...
                st.col.filtered, ...
                'FaceAlpha', 0.10, ...
                'EdgeColor', 'none');

            % Main spectrum on top of the shading.
            hSpec = plot(ax, fPlot, PPlot, ...
                'Color', st.col.noisy, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);

            % Manual dense dash for fc; the built-in '--' looked too sparse.
            dashLength = 0.018 * yRange;
            gapLength = 0.013 * yRange;
            yStart = yLow;
            while yStart < yHigh
                yEnd = min(yStart + dashLength, yHigh);
                plot(ax, [cutoffHz, cutoffHz], [yStart, yEnd], ...
                    'Color', [0, 0, 0], ...
                    'LineStyle', '-', ...
                    'LineWidth', 0.75);
                yStart = yEnd + gapLength;
            end

            % Dummy handle just for the legend label.
            hCut = plot(ax, NaN, NaN, ...
                'Color', [0, 0, 0], ...
                'LineStyle', '--', ...
                'LineWidth', 0.75);

            % Small filled marker for the mains hum peak.
            hHum = plot(ax, humFreq, humMag, 'v', ...
                'LineStyle', 'none', ...
                'MarkerSize', max(st.marker.triangleSize - 1.2, 4.0), ...
                'LineWidth', st.marker.lineWidth, ...
                'MarkerEdgeColor', st.col.peak, ...
                'MarkerFaceColor', st.col.peak);

            % Labels are handled by the legend now.

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, 'DFT with Noise Identification and Cutoff Selection');
            xlabel(ax, 'Frequency (Hz)', st.labelArgs{:});
            ylabel(ax, 'Magnitude (mV)', st.labelArgs{:});

            xlim(ax, [0, 150]);
            set(ax, 'XTick', 0:25:150);
            ylim(ax, [yLow, yHigh]);

            lgd = legend(ax, [hSpec, hHum, hCut, hHF], ...
                {'Noisy ECG spectrum', ...
                 '60 Hz mains hum', ...
                 sprintf('Selected cutoff: %.0f Hz', cutoffHz), ...
                 'High-frequency noise region'}, ...
                'Location', 'northeast', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);
        end


        function fig = filterResponses(F, ~)
            %FILTERRESPONSES Magnitude responses only; phase is not needed for this report.

            st = ecgplot.style();
            nFreq = 4096;

            [HFIR, f] = freqz(F.fir.b, F.fir.a, nFreq, F.fs);
            [HIIR, ~] = freqz(F.iir.b, F.iir.a, nFreq, F.fs);

            magFIR = 20 * log10(max(abs(HFIR), eps));
            magIIR = 20 * log10(max(abs(HIIR), eps));

            fig = ecgplot.newFigure('FIR and IIR Filter Magnitude Responses', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            hFIR = plot(ax, f, magFIR, ...
                'Color', st.col.filtered, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hIIR = plot(ax, f, magIIR, ...
                'Color', st.col.cutoff, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hCut = xline(ax, F.cutoffHz, '--', ...
                sprintf(' f_c = %.0f Hz', F.cutoffHz), ...
                'Color', [0, 0, 0], ...
                'LineWidth', 0.85, ...
                'LabelVerticalAlignment', 'middle', ...
                'LabelHorizontalAlignment', 'left');

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, 'Magnitude Responses of FIR and IIR Low-Pass Filters(5th)');
            xlabel(ax, 'Frequency (Hz)', st.labelArgs{:});
            ylabel(ax, 'Magnitude (dB)', st.labelArgs{:});
            xlim(ax, [0, F.fs/2]);
            ylim(ax, [-80, 5]);
            set(ax, 'XTick', 0:20:(F.fs/2));
            set(hCut, 'FontName', st.font.name, 'FontSize', st.font.legendSize);

            lgd = legend(ax, [hFIR, hIIR], ...
                {'FIR low-pass response', 'IIR Butterworth low-pass response'}, ...
                'Location', 'southwest', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);
        end


        function fig = filteredComparison(t, origSig, noisySig, Y, ~)
            %FILTEREDCOMPARISON Full-length filtered signal comparison.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);
            noisySig = noisySig(:);

            fig = ecgplot.newFigure('Filtered ECG Comparison', 'twoPanel');
            xLimFull = ecgplot.timeLimits(t);
            yLimCommon = ecgplot.calcYLim([origSig; noisySig; Y.fir(:); Y.iir(:)], 0.13);

            ax1 = axes(fig, 'Position', st.axPos.upperPanel);
            hold(ax1, 'on');
            hNoisy1 = plot(ax1, t, noisySig, ...
                'Color', [0.92, 0.67, 0.46], ...
                'LineStyle', '-', ...
                'LineWidth', 0.45);
            hOrig1 = plot(ax1, t, origSig, ...
                'Color', st.col.original, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);
            hFIR = plot(ax1, t, Y.fir, ...
                'Color', st.col.filtered, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);
            hold(ax1, 'off');
            ecgplot.formatAxes(ax1);
            ecgplot.panelTitle(ax1, 'Original ECG Compared with FIR Filter(5th) Output');
            ylabel(ax1, 'Voltage (mV)', st.labelArgs{:});
            xlabel(ax1, 'Time (s)', st.labelArgs{:});
            xlim(ax1, xLimFull);
            ylim(ax1, yLimCommon);
            ecgplot.setSecondTicks(ax1, xLimFull, 1.0);
            lgd1 = legend(ax1, [hOrig1, hFIR, hNoisy1], ...
                {'Original ECG', 'FIR filtered ECG', 'Noisy ECG'}, ...
                'Location', 'southeast', 'Box', 'on', 'Interpreter', 'none');
            ecgplot.formatLegend(lgd1);

            ax2 = axes(fig, 'Position', st.axPos.lowerPanel);
            hold(ax2, 'on');
            hNoisy2 = plot(ax2, t, noisySig, ...
                'Color', [0.92, 0.67, 0.46], ...
                'LineStyle', '-', ...
                'LineWidth', 0.45);
            hOrig2 = plot(ax2, t, origSig, ...
                'Color', st.col.original, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.signal);
            hIIR = plot(ax2, t, Y.iir, ...
                'Color', st.col.cutoff, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);
            hold(ax2, 'off');
            ecgplot.formatAxes(ax2);
            ecgplot.panelTitle(ax2, 'Original ECG Compared with IIR Butterworth Filter(5th) Output');
            ylabel(ax2, 'Voltage (mV)', st.labelArgs{:});
            xlabel(ax2, 'Time (s)', st.labelArgs{:});
            xlim(ax2, xLimFull);
            ylim(ax2, yLimCommon);
            ecgplot.setSecondTicks(ax2, xLimFull, 1.0);
            lgd2 = legend(ax2, [hOrig2, hIIR, hNoisy2], ...
                {'Original ECG', 'IIR filtered ECG', 'Noisy ECG'}, ...
                'Location', 'southeast', 'Box', 'on', 'Interpreter', 'none');
            ecgplot.formatLegend(lgd2);

            linkaxes([ax1, ax2], 'x');
        end


        function fig = filteredComparisonZoom(t, origSig, noisySig, Y, zoomRange, ~)
            %FILTEREDCOMPARISONZOOM Local filtered comparison around one clearer beat.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);
            noisySig = noisySig(:);
            idx = t >= zoomRange(1) & t <= zoomRange(2);

            fig = ecgplot.newFigure('Zoomed Filtered ECG Comparison', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            noisyCol = ecgplot.blendWithWhite([0.92, 0.67, 0.46], 0.45);
            origCol  = ecgplot.blendWithWhite(st.col.original, 0.82);
            firCol   = ecgplot.blendWithWhite(st.col.filtered, 0.82);
            iirCol   = ecgplot.blendWithWhite(st.col.cutoff, 0.82);

            hNoisy = plot(ax, t(idx), noisySig(idx), ...
                'Color', noisyCol, ...
                'LineStyle', '-', ...
                'LineWidth', 0.55);

            hFIR = plot(ax, t(idx), Y.fir(idx), ...
                'Color', firCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hIIR = plot(ax, t(idx), Y.iir(idx), ...
                'Color', iirCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hOrig = plot(ax, t(idx), origSig(idx), ...
                'Color', origCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, sprintf('Filtered ECG Detail, %.2f--%.2f s', zoomRange(1), zoomRange(2)));
            xlabel(ax, 'Time (s)', st.labelArgs{:});
            ylabel(ax, 'Voltage (mV)', st.labelArgs{:});
            xlim(ax, zoomRange);
            yMain = ecgplot.calcYLim([origSig(idx); noisySig(idx); Y.fir(idx); Y.iir(idx)], 0.16);
            ylim(ax, yMain);
            ecgplot.setSecondTicks(ax, zoomRange, 0.10);

            lgd = legend(ax, [hOrig, hFIR, hIIR, hNoisy], ...
                {'Original ECG', 'FIR filtered ECG', 'IIR filtered ECG', 'Noisy ECG'}, ...
                'Location', 'southeast', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);
        end


        function fig = filteredSpectrumComparison(noisySig, Y, fs, cutoffHz, ~)
            %FILTEREDSPECTRUMCOMPARISON Frequency-domain check after filtering.
            %
            % The main axis keeps the 0--150 Hz overview; the inset is just to make
            % the 60 Hz reduction visible without zooming the whole figure.

            st = ecgplot.style();

            [f, PNoisy] = ecgplot.singleSidedSpectrum(noisySig, fs);
            [~, PFIR]   = ecgplot.singleSidedSpectrum(Y.fir, fs);
            [~, PIIR]   = ecgplot.singleSidedSpectrum(Y.iir, fs);

            plotMask = f <= 150;
            fPlot = f(plotMask);
            PNoisyPlot = PNoisy(plotMask);
            PFIRPlot = PFIR(plotMask);
            PIIRPlot = PIIR(plotMask);

            yMax = max([PNoisyPlot(:); PFIRPlot(:); PIIRPlot(:)]);
            if ~isfinite(yMax) || yMax <= 0
                yMax = 1;
            end

            yPadTop = 0.13 * yMax;
            yPadBottom = 0.03 * yMax;
            yLow = -yPadBottom;
            yHigh = yMax + yPadTop;
            yRange = yHigh - yLow;

            % Re-find the noisy 60 Hz peak so the inset sits in the right place.
            humMask = f >= 57 & f <= 63;
            humFreq = 60;
            humMag = interp1(f, PNoisy, humFreq, 'linear', 'extrap');
            if any(humMask)
                fHum = f(humMask);
                PHum = PNoisy(humMask);
                [humMag, idxHum] = max(PHum);
                humFreq = fHum(idxHum);
            end

            % Same HF artificial-noise band as the annotated DFT.
            hfBand = [117, 142];

            % I softened the colours by blending with white; line alpha is messy across MATLAB versions.
            noisyCol = ecgplot.blendWithWhite([0.92, 0.67, 0.46], 0.55);
            firCol   = ecgplot.blendWithWhite(st.col.filtered, 0.62);
            iirCol   = ecgplot.blendWithWhite(st.col.cutoff,   0.58);

            fig = ecgplot.newFigure('Filtered Frequency-Domain Comparison', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            % Shade first so it stays behind the spectrum.
            hHF = patch(ax, ...
                [hfBand(1), hfBand(2), hfBand(2), hfBand(1)], ...
                [yLow, yLow, yHigh, yHigh], ...
                st.col.filtered, ...
                'FaceAlpha', 0.10, ...
                'EdgeColor', 'none');

            hNoisy = plot(ax, fPlot, PNoisyPlot, ...
                'Color', noisyCol, ...
                'LineStyle', '-', ...
                'LineWidth', 0.80);

            hFIR = plot(ax, fPlot, PFIRPlot, ...
                'Color', firCol, ...
                'LineStyle', '-', ...
                'LineWidth', 1.05);

            hIIR = plot(ax, fPlot, PIIRPlot, ...
                'Color', iirCol, ...
                'LineStyle', '-', ...
                'LineWidth', 1.00);

            % Dense fc line again, to match the DFT figure.
            dashLength = 0.018 * yRange;
            gapLength = 0.013 * yRange;
            yStart = yLow;
            while yStart < yHigh
                yEnd = min(yStart + dashLength, yHigh);
                plot(ax, [cutoffHz, cutoffHz], [yStart, yEnd], ...
                    'Color', [0, 0, 0], ...
                    'LineStyle', '-', ...
                    'LineWidth', 0.75);
                yStart = yEnd + gapLength;
            end

            hCut = plot(ax, NaN, NaN, ...
                'Color', [0, 0, 0], ...
                'LineStyle', '--', ...
                'LineWidth', 0.75);

            hHum = plot(ax, humFreq, humMag, 'v', ...
                'LineStyle', 'none', ...
                'MarkerSize', max(st.marker.triangleSize - 1.2, 4.0), ...
                'LineWidth', st.marker.lineWidth, ...
                'MarkerEdgeColor', st.col.peak, ...
                'MarkerFaceColor', st.col.peak);

            % Box shows where the inset is coming from.
            zoomBand = [57, 63];
            zoomMaskMain = fPlot >= zoomBand(1) & fPlot <= zoomBand(2);
            zoomValsMain = [PNoisyPlot(zoomMaskMain); PFIRPlot(zoomMaskMain); PIIRPlot(zoomMaskMain)];
            zLow = min(zoomValsMain);
            zHigh = max(zoomValsMain);
            zRange = zHigh - zLow;
            if ~isfinite(zRange) || zRange <= 0
                zRange = max(abs(zHigh), 1);
            end
            zoomYLow = max(yLow, zLow - 0.25 * zRange);
            zoomYHigh = min(yHigh, zHigh + 0.35 * zRange);
            rectangle(ax, ...
                'Position', [zoomBand(1), zoomYLow, diff(zoomBand), zoomYHigh - zoomYLow], ...
                'EdgeColor', [0.22, 0.22, 0.22], ...
                'LineStyle', '-', ...
                'LineWidth', 0.75);

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, 'Frequency-Domain Comparison');
            xlabel(ax, 'Frequency (Hz)', st.labelArgs{:});
            ylabel(ax, 'Magnitude (mV)', st.labelArgs{:});
            xlim(ax, [0, 150]);
            set(ax, 'XTick', 0:25:150);
            ylim(ax, [yLow, yHigh]);

            % Legend kept in the corner to avoid covering the 60 Hz inset.
            lgd = legend(ax, [hNoisy, hFIR, hIIR, hHum, hCut, hHF], ...
                {'Noisy ECG spectrum', ...
                 'FIR filtered spectrum', ...
                 'IIR filtered spectrum', ...
                 '60 Hz mains hum', ...
                 sprintf('Selected cutoff: %.0f Hz', cutoffHz), ...
                 'High-frequency noise region'}, ...
                'Location', 'northeast', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);

            % Compact inset: big enough to show 60 Hz, but not large enough to dominate.
            axInset = axes(fig, 'Position', [0.41, 0.4, 0.25, 0.23]);
            hold(axInset, 'on');

            insetMask = f >= zoomBand(1) & f <= zoomBand(2);
            plot(axInset, f(insetMask), PNoisy(insetMask), ...
                'Color', noisyCol, ...
                'LineStyle', '-', ...
                'LineWidth', 0.80);
            plot(axInset, f(insetMask), PFIR(insetMask), ...
                'Color', firCol, ...
                'LineStyle', '-', ...
                'LineWidth', 1.00);
            plot(axInset, f(insetMask), PIIR(insetMask), ...
                'Color', iirCol, ...
                'LineStyle', '-', ...
                'LineWidth', 0.95);

            plot(axInset, humFreq, humMag, 'v', ...
                'LineStyle', 'none', ...
                'MarkerSize', max(st.marker.triangleSize - 2.0, 3.8), ...
                'LineWidth', st.marker.lineWidth, ...
                'MarkerEdgeColor', st.col.peak, ...
                'MarkerFaceColor', st.col.peak);

            plot(axInset, [60, 60], [0, max([PNoisy(insetMask); PFIR(insetMask); PIIR(insetMask)])], ...
                'Color', [0.22, 0.22, 0.22], ...
                'LineStyle', ':', ...
                'LineWidth', 0.70);

            hold(axInset, 'off');
            ecgplot.formatAxes(axInset);
            title(axInset, 'Zoom at 60 Hz', ...
                'FontName', st.font.name, ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                'Interpreter', 'none');
            xlim(axInset, zoomBand);
            insetVals = [PNoisy(insetMask); PFIR(insetMask); PIIR(insetMask)];
            iLow = min(insetVals);
            iHigh = max(insetVals);
            iRange = iHigh - iLow;
            if ~isfinite(iRange) || iRange <= 0
                iRange = max(abs(iHigh), 1);
            end
            ylim(axInset, [iLow - 0.20 * iRange, iHigh + 0.35 * iRange]);
            set(axInset, ...
                'FontSize', 7.5, ...
                'XTick', 55:5:65, ...
                'YTickMode', 'auto', ...
                'Box', 'on', ...
                'Layer', 'top');
            xlabel(axInset, 'Hz', 'FontName', st.font.name, 'FontSize', 7.5);
            ylabel(axInset, 'mV', 'FontName', st.font.name, 'FontSize', 7.5);
        end


        function fig = improvementComparison(t, origSig, noisySig, Y, zoomRange, cfg)
            %IMPROVEMENTCOMPARISON Standard IIR versus the extra cleaned-up signal.
            %
            % This is the Step 11 figure: notch the 60 Hz line, then use zero-phase
            % wider LPF so small P/T features are less smeared.

            st = ecgplot.style();
            t = t(:);
            origSig = origSig(:);
            noisySig = noisySig(:);

            if ~isfield(Y, 'improved')
                error('Y.improved is missing. Run ecgtool.applyFilters with the improved filter design.');
            end

            idx = t >= zoomRange(1) & t <= zoomRange(2);

            fig = ecgplot.newFigure('Improved ECG Filter Comparison', 'wide');
            ax = axes(fig, 'Position', st.axPos.single);
            hold(ax, 'on');

            noisyCol    = ecgplot.blendWithWhite([0.92, 0.67, 0.46], 0.45);
            origCol     = ecgplot.blendWithWhite(st.col.original, 0.92);
            improvedCol = ecgplot.blendWithWhite(st.col.improved, 0.64);
            iirCol      = ecgplot.blendWithWhite(st.col.cutoff, 0.64);

            hNoisy = plot(ax, t(idx), noisySig(idx), ...
                'Color', noisyCol, ...
                'LineStyle', '-', ...
                'LineWidth', 0.55);

            hImproved = plot(ax, t(idx), Y.improved(idx), ...
                'Color', improvedCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            hIIR = plot(ax, t(idx), Y.iir(idx), ...
                'Color', iirCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);

            % Draw original later; otherwise the comparison trace gets buried.
            hOrig = [];

            hold(ax, 'off');

            ecgplot.formatAxes(ax);
            ecgplot.panelTitle(ax, sprintf('Improvement Performance, %.1f--%.1f s', zoomRange(1), zoomRange(2)));
            xlabel(ax, 'Time (s)', st.labelArgs{:});
            ylabel(ax, 'Voltage (mV)', st.labelArgs{:});
            xlim(ax, zoomRange);
            yMain = ecgplot.calcYLim([origSig(idx); noisySig(idx); Y.iir(idx); Y.improved(idx)], 0.16);
            ylim(ax, yMain);
            ecgplot.setSecondTicks(ax, zoomRange, 0.10);

            % Small P/T insets, because these parts disappear in the full-scale view.
            tWin = t(idx);
            xOrigMain = origSig(idx);
            [~, iR] = max(xOrigMain);
            rTime = tWin(iR);

            pRange = [max(zoomRange(1), rTime - 0.200), min(zoomRange(2), rTime - 0.070)];
            tRange = [max(zoomRange(1), rTime + 0.170), min(zoomRange(2), rTime + 0.460)];

            % Mark where the insets come from on the main axis as well.
            pMaskMain = t >= pRange(1) & t <= pRange(2);
            tMaskMain = t >= tRange(1) & t <= tRange(2);

            yP = [origSig(pMaskMain); Y.iir(pMaskMain); Y.improved(pMaskMain)];
            yT = [origSig(tMaskMain); Y.iir(tMaskMain); Y.improved(tMaskMain)];

            pYL = ecgplot.calcYLim(yP, 0.12);
            tYL = ecgplot.calcYLim(yT, 0.12);

            rectangle(ax, ...
                'Position', [pRange(1), pYL(1), diff(pRange), diff(pYL)], ...
                'EdgeColor', [0.20, 0.20, 0.20], ...
                'LineStyle', '-', ...
                'LineWidth', 0.95, ...
                'FaceColor', 'none');

            rectangle(ax, ...
                'Position', [tRange(1), tYL(1), diff(tRange), diff(tYL)], ...
                'EdgeColor', [0.20, 0.20, 0.20], ...
                'LineStyle', '-', ...
                'LineWidth', 0.95, ...
                'FaceColor', 'none');

            text(ax, pRange(1), pYL(2), ' P', ...
                'FontName', st.font.name, ...
                'FontSize', st.font.legendSize, ...
                'FontWeight', 'bold', ...
                'Color', [0.15, 0.15, 0.15], ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');

            text(ax, tRange(1), tYL(2), ' T', ...
                'FontName', st.font.name, ...
                'FontSize', st.font.legendSize, ...
                'FontWeight', 'bold', ...
                'Color', [0.15, 0.15, 0.15], ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');

            % Blue reference goes last so it stays visible.
            hold(ax, 'on');
            hOrig = plot(ax, t(idx), origSig(idx), ...
                'Color', origCol, ...
                'LineStyle', '-', ...
                'LineWidth', st.line.emphasis);
            hold(ax, 'off');

            axInsetP = axes(fig, 'Position', [0.58, 0.66, 0.30, 0.15]);
            hold(axInsetP, 'on');
            pMask = t >= pRange(1) & t <= pRange(2);
            plot(axInsetP, t(pMask), Y.improved(pMask), 'Color', improvedCol, 'LineStyle', '-', 'LineWidth', 1.05);
            plot(axInsetP, t(pMask), Y.iir(pMask), 'Color', iirCol, 'LineStyle', '-', 'LineWidth', 1.05);
            plot(axInsetP, t(pMask), origSig(pMask), 'Color', origCol, 'LineStyle', '-', 'LineWidth', 1.15);
            hold(axInsetP, 'off');
            ecgplot.formatAxes(axInsetP);
            title(axInsetP, 'P-wave detail', 'FontName', st.font.name, 'FontSize', st.font.legendSize, 'FontWeight', 'normal');
            xlim(axInsetP, pRange);
            ylim(axInsetP, ecgplot.calcYLim([origSig(pMask); Y.iir(pMask); Y.improved(pMask)], 0.14));
            set(axInsetP, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on', 'Color', 'white', 'Layer', 'top');

            axInsetT = axes(fig, 'Position', [0.58, 0.46, 0.30, 0.15]);
            hold(axInsetT, 'on');
            tMask2 = t >= tRange(1) & t <= tRange(2);
            plot(axInsetT, t(tMask2), Y.improved(tMask2), 'Color', improvedCol, 'LineStyle', '-', 'LineWidth', 1.05);
            plot(axInsetT, t(tMask2), Y.iir(tMask2), 'Color', iirCol, 'LineStyle', '-', 'LineWidth', 1.05);
            plot(axInsetT, t(tMask2), origSig(tMask2), 'Color', origCol, 'LineStyle', '-', 'LineWidth', 1.15);
            hold(axInsetT, 'off');
            ecgplot.formatAxes(axInsetT);
            title(axInsetT, 'T-wave detail', 'FontName', st.font.name, 'FontSize', st.font.legendSize, 'FontWeight', 'normal');
            xlim(axInsetT, tRange);
            ylim(axInsetT, ecgplot.calcYLim([origSig(tMask2); Y.iir(tMask2); Y.improved(tMask2)], 0.14));
            set(axInsetT, 'XTickLabel', [], 'YTickLabel', [], 'Box', 'on', 'Color', 'white', 'Layer', 'top');

            % Do not call axes(ax) here, otherwise MATLAB can cover the insets.

            if isstruct(cfg) && isfield(cfg, 'notchFreqHz') && isfield(cfg, 'improvedCutoffHz')
                improvedLabel = sprintf('Improved: %.0f Hz notch + %.0f Hz zero-phase LP', ...
                    cfg.notchFreqHz, cfg.improvedCutoffHz);
            else
                improvedLabel = 'Improved notch + zero-phase LP';
            end

            lgd = legend(ax, [hOrig, hIIR, hImproved, hNoisy], ...
                {'Original ECG', 'Standard IIR filtered ECG', improvedLabel, 'Noisy ECG'}, ...
                'Location', 'southeast', ...
                'Box', 'on', ...
                'Interpreter', 'none');
            ecgplot.formatLegend(lgd);

            % Keep the inset axes above the main plot after the legend is made.
            try
                uistack(axInsetP, 'top');
                uistack(axInsetT, 'top');
            catch
                % Some MATLAB installs do not expose uistack; the plot is still OK without it.
            end
        end


        function fig = performanceMetricsBars(metricsTable)
            %PERFORMANCEMETRICSBARS Compact bar chart for the key waveform metrics.
            %
            % I left the spectral/SNR metrics in the CSV table and kept this figure
            % focused on what is easiest to explain visually: error, delay and PQRST
            % shape preservation.

            st = ecgplot.style();

            % Keep the three filter rows only, in the order used in the report.
            wanted = {'FIR', 'IIR', 'Improved'};
            rowIdx = zeros(1, numel(wanted));

            filterClass = metricsTable.Filter_Class;
            if iscategorical(filterClass)
                filterClass = cellstr(filterClass);
            elseif isstring(filterClass)
                filterClass = cellstr(filterClass);
            end

            for k = 1:numel(wanted)
                idx = find(strcmp(filterClass, wanted{k}), 1, 'first');
                if isempty(idx)
                    error('Could not find row for Filter_Class = %s.', wanted{k});
                end
                rowIdx(k) = idx;
            end

            T = metricsTable(rowIdx, :);

            labels = {'FIR', 'IIR', 'Improved'};
            barColors = [ ...
                st.col.filtered; ...
                st.col.cutoff; ...
                st.col.improved ...
            ];

            metricVals = { ...
                T.RMSE_mV, ...
                T.R_abs_delay_ms, ...
                T.PQRST_RMSE_mV, ...
                T.PQRST_Corr ...
            };

            metricTitles = { ...
                'Overall RMSE', ...
                'R-peak Delay', ...
                'PQRST RMSE', ...
                'PQRST Correlation' ...
            };

            yLabels = { ...
                'mV', ...
                'ms', ...
                'mV', ...
                'Correlation' ...
            };

            directionText = { ...
                '', ...
                '', ...
                '', ...
                '' ...
            };

            % Short wide layout so it fits under the text more easily.
            fig = ecgplot.newFigure('Key PQRST Metrics Comparison', 'wide');
            set(fig, ...
                'Position', [1.0, 1.0, 8.35, 2.85], ...
                'PaperUnits', 'inches', ...
                'PaperPositionMode', 'manual', ...
                'PaperPosition', [0, 0, 8.35, 2.85], ...
                'PaperSize', [8.35, 2.85], ...
                'InvertHardcopy', 'off', ...
                'Renderer', 'painters');

            % Manual 1x4 positions export more predictably than subplot here.
            axPos = [ ...
                0.060, 0.245, 0.205, 0.585; ...
                0.300, 0.245, 0.205, 0.585; ...
                0.540, 0.245, 0.205, 0.585; ...
                0.780, 0.245, 0.205, 0.585 ...
            ];

            for k = 1:4
                ax = axes(fig, 'Position', axPos(k, :));

                vals = double(metricVals{k}(:)).';
                b = bar(ax, vals, 0.68, 'FaceColor', 'flat');
                b.CData = barColors;

                ecgplot.formatAxes(ax);
                title(ax, metricTitles{k}, ...
                    'FontName', st.font.name, ...
                    'FontSize', 10, ...
                    'FontWeight', 'bold', ...
                    'Interpreter', 'none');

                ylabel(ax, yLabels{k}, st.labelArgs{:});
                set(ax, 'XTick', 1:3, 'XTickLabel', labels);
                xtickangle(ax, 0);

                yMin = min([0, vals]);
                yMax = max([0, vals]);
                ySpan = yMax - yMin;
                if ySpan <= 0
                    ySpan = 1;
                end

                if strcmp(metricTitles{k}, 'PQRST Correlation')
                    % Correlation has a hard upper bound, so leave room for the label.
                    ylim(ax, [0, 1.18]);
                    ySpanForText = 1;
                else
                    % More top padding because the numeric labels sit above bars.
                    padBottom = 0.12 * ySpan;
                    padTop = 0.32 * ySpan;
                    ylim(ax, [yMin - padBottom, yMax + padTop]);
                    ySpanForText = ySpan;
                end

                % Small panel note instead of an overlong title.
                text(ax, 0.04, 0.90, directionText{k}, ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top', ...
                    'FontName', st.font.name, ...
                    'FontSize', 8, ...
                    'Color', [0.25, 0.25, 0.25]);

                % Numbers on the bars make the figure easier to quote in the report.
                for n = 1:numel(vals)
                    if vals(n) >= 0
                        yTxt = vals(n) + 0.040 * ySpanForText;
                        vAlign = 'bottom';
                    else
                        yTxt = vals(n) - 0.040 * ySpanForText;
                        vAlign = 'top';
                    end

                    if strcmp(metricTitles{k}, 'PQRST Correlation')
                        txt = sprintf('%.4f', vals(n));
                    else
                        txt = sprintf('%.2f', vals(n));
                    end

                    text(ax, n, yTxt, txt, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', vAlign, ...
                        'FontName', st.font.name, ...
                        'FontSize', 8, ...
                        'Color', st.col.text);
                end
            end
        end


        function saveFigure(fig, baseFile)
            %SAVEFIGURE Save one figure as PNG and PDF.
            % I keep this as a tiny wrapper so main.m only passes 'Fig1', 'Fig2', etc.

            [outDir, ~, ~] = fileparts(baseFile);
            if ~isempty(outDir) && ~exist(outDir, 'dir')
                mkdir(outDir);
            end

            pngFile = [baseFile, '.png'];
            pdfFile = [baseFile, '.pdf'];

            try
                exportgraphics(fig, pngFile, 'Resolution', 300);
                exportgraphics(fig, pdfFile, 'ContentType', 'vector');
            catch
                % Fallback for older MATLAB versions.
                set(fig, 'PaperPositionMode', 'auto');
                print(fig, pngFile, '-dpng', '-r300');
                print(fig, pdfFile, '-dpdf', '-painters');
            end
        end


        function fig = newFigure(figName, layoutType)
            %NEWFIGURE Create a clean white figure with the chosen size.

            st = ecgplot.style();
            fig = figure('Color', 'w', 'Name', figName, 'Units', 'inches');

            switch lower(layoutType)
                case 'twopanel'
                    pos = [1.0, 1.0, st.fig.doubleColumnWidth, st.fig.twoPanelHeight];
                case 'wide'
                    pos = [1.0, 1.0, st.fig.doubleColumnWidth, st.fig.wideHeight];
                otherwise
                    pos = [1.0, 1.0, st.fig.singleColumnWidth, st.fig.singleHeight];
            end

            set(fig, ...
                'Position', pos, ...
                'PaperUnits', 'inches', ...
                'PaperPositionMode', 'manual', ...
                'PaperPosition', [0, 0, pos(3), pos(4)], ...
                'PaperSize', [pos(3), pos(4)], ...
                'InvertHardcopy', 'off', ...
                'Renderer', 'painters');
        end


        function formatAxes(ax)
            %FORMATAXES My default axes formatting.

            st = ecgplot.style();

            set(ax, ...
                'FontName', st.font.name, ...
                'FontSize', st.font.tickSize, ...
                'TickLabelInterpreter', 'none', ...
                'LineWidth', st.axes.lineWidth, ...
                'TickDir', 'out', ...
                'TickLength', st.axes.tickLength, ...
                'Box', 'on', ...
                'Layer', 'top', ...
                'XMinorTick', 'on', ...
                'YMinorTick', 'on', ...
                'XGrid', 'on', ...
                'YGrid', 'on', ...
                'XMinorGrid', 'off', ...
                'YMinorGrid', 'off', ...
                'GridLineStyle', '-', ...
                'GridColor', st.col.grid, ...
                'GridAlpha', st.axes.gridAlpha);
        end


        function panelTitle(ax, titleText)
            %PANELTITLE Standard title formatting.

            st = ecgplot.style();

            title(ax, titleText, ...
                'FontName', st.font.name, ...
                'FontSize', st.font.titleSize, ...
                'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end


        function formatLegend(lgd)
            %FORMATLEGEND Standard legend formatting.

            st = ecgplot.style();
            set(lgd, ...
                'FontName', st.font.name, ...
                'FontSize', st.font.legendSize, ...
                'TextColor', st.col.text, ...
                'Color', 'white', ...
                'EdgeColor', st.col.legendEdge, ...
                'LineWidth', st.legend.lineWidth);
        end


        function yLim = calcYLim(yData, frac)
            %CALCYLIM Add a bit of y padding so traces do not touch the box.

            yData = double(yData(:));
            yData = yData(isfinite(yData));

            if isempty(yData)
                yLim = [0, 1];
                return;
            end

            ymin = min(yData);
            ymax = max(yData);
            yrange = ymax - ymin;

            if yrange == 0
                yrange = max(abs(ymax), 1);
            end

            yLim = [ymin - frac * yrange, ymax + frac * yrange];
        end


        function padY(ax, yData, frac)
            %PADY Old wrapper kept so earlier plot calls still work.

            ylim(ax, ecgplot.calcYLim(yData, frac));
        end


        function xLim = timeLimits(t)
            %TIMELIMITS Prefer clean second boundaries.

            t = double(t(:));
            t = t(isfinite(t));

            if isempty(t)
                xLim = [0, 1];
                return;
            end

            xmin = min(t);
            xmax = max(t);

            if xmin >= 0 && xmax > 1
                xmin = 0;
                xmax = ceil(xmax);
            end

            xLim = [xmin, xmax];
        end


        function setSecondTicks(ax, xLim, step)
            %SETSECONDTICKS Clean time ticks if they will not overcrowd.

            if step <= 0
                return;
            end

            ticks = xLim(1):step:xLim(2);
            if numel(ticks) >= 2 && numel(ticks) <= 16
                set(ax, 'XTick', ticks);
            end
        end



        function hLegend = denseDashedPlot(ax, x, y, lineColor, lineWidth, dashLen, gapLen)
            %DENSEDASHEDPLOT Draw dense dashes manually.
            % MATLAB does not give nice control over '--' spacing, so I segment it.

            x = x(:);
            y = y(:);
            valid = isfinite(x) & isfinite(y);
            x = x(valid);
            y = y(valid);

            if nargin < 6 || isempty(dashLen)
                dashLen = 4;
            end
            if nargin < 7 || isempty(gapLen)
                gapLen = 2;
            end

            wasHold = ishold(ax);
            hold(ax, 'on');
            k = 1;
            while k <= numel(x)
                k2 = min(k + dashLen - 1, numel(x));
                plot(ax, x(k:k2), y(k:k2), ...
                    'Color', lineColor, ...
                    'LineStyle', '-', ...
                    'LineWidth', lineWidth);
                k = k2 + gapLen + 1;
            end

            hLegend = plot(ax, NaN, NaN, ...
                'Color', lineColor, ...
                'LineStyle', '--', ...
                'LineWidth', lineWidth);

            if ~wasHold
                hold(ax, 'off');
            end
        end


        function cOut = blendWithWhite(cIn, alphaLike)
            %BLENDWITHWHITE Fake transparency by mixing a colour with white.
            % alphaLike = 1 keeps the colour; alphaLike = 0 gives white.

            cIn = double(cIn(:)).';
            alphaLike = max(0, min(1, double(alphaLike)));
            cOut = alphaLike * cIn + (1 - alphaLike) * [1, 1, 1];
        end

        function st = style()
            %STYLE Central plotting style used by every figure.
            %
            % Most visual edits are here.  The margins are slightly generous on
            % purpose because exported PDF labels were getting clipped before.

            st = struct();

            st.fig.singleColumnWidth = 3.50;
            st.fig.doubleColumnWidth = 8.35;
            st.fig.singleHeight = 3.20;
            st.fig.wideHeight = 4.95;
            st.fig.twoPanelHeight = 7.40;

            st.font.name = 'Times New Roman';
            st.font.tickSize = 9;
            st.font.labelSize = 10;
            st.font.titleSize = 11;
            st.font.legendSize = 9;

            st.axes.lineWidth = 0.75;
            st.axes.tickLength = [0.007, 0.007];
            st.axes.gridAlpha = 0.18;

            st.legend.lineWidth = 0.75;

            % These margins are a little conservative, but they stop PDF clipping.
            st.axPos.single = [0.118, 0.165, 0.824, 0.710];
            st.axPos.upperPanel = [0.118, 0.610, 0.824, 0.285];
            st.axPos.lowerPanel = [0.118, 0.160, 0.824, 0.285];

            st.line.noisy = 0.70;
            st.line.signal = 1.10;
            st.line.emphasis = 1.45;
            st.line.zoomOriginal = 1.18;

            st.marker.size = 8.0;
            st.marker.triangleSize = 6.2;
            st.marker.lineWidth = 1.20;

            % Palette based on Okabe-Ito; still adjusted a bit for my overlays.
            st.col.original   = [0.000, 0.260, 0.520];   % navy blue
            st.col.originalSoft = [0.180, 0.430, 0.680]; % lighter blue for overlays
            st.col.noisy      = [0.835, 0.369, 0.000];   % muted vermillion
            st.col.filtered   = [0.000, 0.500, 0.500];   % teal
            st.col.improved   = [0.000, 0.410, 0.160];   % dark green for Step 11 improvement
            st.col.peak       = [0.700, 0.070, 0.170];   % crimson
            st.col.cutoff     = [0.494, 0.184, 0.556];   % purple
            st.col.grid       = [0.865, 0.865, 0.865];
            st.col.legendEdge = [0.250, 0.250, 0.250];
            st.col.text       = [0.000, 0.000, 0.000];

            st.labelArgs = { ...
                'FontName', st.font.name, ...
                'FontSize', st.font.labelSize, ...
                'FontWeight', 'normal', ...
                'Interpreter', 'none'};
        end

    end
end
