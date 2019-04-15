% Calculate fft for binary sequences to get strongest frequencies
% (correspoding to rhythmicity)

% source: http://de.mathworks.com/help/matlab/ref/fft.html
function FFT
Fs = 200;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 861;             % Length of signal
t = (0:L-1)*T;        % Time vector
binarydata = binarydata
for i = 1:length(binarydata)
X = binarydata(:,i);

Y = fft(X);

% for i = 1: length (Y)
%     Y (i) = sum (Y (i,:));
% end
% Y = Y';
%Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Fs/L;
%f = [ans:ans*2:100.2326]; %0-100 Hz

f = Fs*(0:(L/2))/L; % 0-50 Hz, Original

% OUTPUT: List of strongest responses with corresponding frequency
% (descending order to easily get "best" pulse frequency
Optimum (:,1) = P1;
Optimum (:,2) = f;
optimal_frequencies = sortrows (Optimum, 1);

%Delet all values for f < 5
 indices = find(abs(Optimum(:,2)<5));
    Optimum(indices,:) = NaN; 
 %find f for best P1
[maxP1,ind] = max(Optimum(:,1))
bestf(i) = Optimum(ind,2)

end 
%Plot Frequenz Domain
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% xlim([0 100])
end