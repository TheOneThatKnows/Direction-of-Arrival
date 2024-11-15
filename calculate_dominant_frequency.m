function [dominant_freq] = calculate_dominant_frequency(signal, fs)
    % Input parameters:
    % signal: Input sequence/signal
    % fs: Sampling frequency in Hz
    %
    % Output:
    % dominant_freq: Dominant frequency in Hz
    
    % Remove DC component (mean)
    signal = signal - mean(signal);
    
    % Calculate FFT
    N = length(signal);
    Y = fft(signal);
    
    % Take only the first half (due to symmetry of FFT)
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Create frequency array
    f = fs*(0:(N/2))/N;
    
    % Find the peak frequency
    [~, idx] = max(P1);
    dominant_freq = f(idx);
    
    % Optional: Plot the frequency spectrum
    figure;
    plot(f, P1);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Single-Sided Amplitude Spectrum');
    hold on;
    plot(dominant_freq, P1(idx), 'r*', 'MarkerSize', 10);
    text(dominant_freq, P1(idx), ['\leftarrow ', num2str(dominant_freq), ' Hz'], ...
        'VerticalAlignment', 'bottom');
end
