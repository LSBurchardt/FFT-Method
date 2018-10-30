%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFT ANALYSIS
%author: Lara S. Burchardt, l.s.burchardt@gmx.de
%version: 0.1.0 
%12.09.2018
%combined originally from: binary_2,FFT_3,FFT_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fft_analysis_master %main function

SetPath                     %subfunction 1, labeled at beginning and end
binary                      %subfunction 2, labeled at beginning and end
FFT                         %subfunction 3, labeled at beginning and end
%gof_original                %subfunction 4, labeled at beginning and end
%random_sequence             %subfunction 5, labeled at beginning and end
%fft_model                   %subfunction 6, labeled at beginning and end
%gof_model                   %subfunction 7, labeled at beginning and end

%% 01: GLOBALS
% definde global variables

FS = 200;               % Sampling frequency                    
T = 1/FS;               % Sampling period       
%t = (0:L-1)*T;          % Time vector

%% 02: DATA MANAGMENT
% load data
% Initialization steps:
%clc;    % Clear the command window.
function SetPath        %subfunction 1

workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;

% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox');
if ~exist(start_path, 'dir')
	start_path = matlabroot;
end
% Ask user to confirm the folder, or change it.
uiwait(msgbox('Pick a starting folder on the next window that will come up.'));
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
fprintf('The top level folder is "%s".\n', topLevelFolder);

% Specify the file pattern.
% Get ALL files using the pattern *.*
% Note the special file pattern.  It has /**/ in it if you want to get files in subfolders of the top level folder.
% filePattern = sprintf('%s/**/*.m; %s/**/*.xml', topLevelFolder, topLevelFolder);
filePattern = sprintf('%s/**/*.xlsx', topLevelFolder); %for .xlsx data
%filePattern = sprintf('%s/**/*.xls', topLevelFolder); %for .xls data
allFileInfo = dir(filePattern);

% Throw out any folders.  We want files only, not folders.
isFolder = [allFileInfo.isdir]; % Logical list of what item is a folder or not.
% Now set those folder entries to null, essentially deleting/removing them from the list.
allFileInfo(isFolder) = [];
% Get a cell array of strings.  
listOfFolderNames = unique({allFileInfo.folder});
numberOfFolders = length(listOfFolderNames);
fprintf('The total number of folders to look in is %d.\n', numberOfFolders);

% Step two: Get a cell array of base filename strings. 
listOfFileNames = {};
listOfFileNames = {allFileInfo.name};
totalNumberOfFiles = length(listOfFileNames);
fprintf('The total number of files in those %d folders is %d.\n', numberOfFolders, totalNumberOfFiles);
end                         % subfunction 1


%% 03: TRANSFORM TO BINARY
% timestamps are converted into a binary sequence

function binary             % subfunction 2

%Step three: loop through files, read .xls files (xlsread) and calculate binary
%code for syllable onsets

binarydata = []
binarydata(:,length(listOfFileNames)) = 0;

for k= 1: length(listOfFileNames)
    

    matfilename = listOfFileNames{:,k};
   
    data = xlsread(matfilename); 
    vectorA = data (:,1); % we just need the syllable onsets 
    vectorA = round(vectorA * 1000 / 5); %we multiply by 1000 to get data 
    % in miliseconds, divide by 5/1 because that is the timeresolution i want
    % to have and round to integer (200 Hz sampling rate--> sample every 5 milliseconds, sampling rate 1000 Hz--> sample every milisecond)

%     binarydata = [1: max(vectorA)];
%     binarydata = binarydata';
%     binarydata(:) = 0;
    
    for l = 1:max(vectorA)
        for i = 1:length(vectorA)
            if l == vectorA(i, 1)
            binarydata(l,k) = 1 ; 
            %else binarydata (l,k) = 0;
            end 
        end 
    end
    
 
end

xlswrite('S_bil_pups_binary_fs200.xlsx', binarydata )
save('S_bil_pups_binary_fs200.mat')

end                             % subfunction 2



%% 04.1: FFT
% source: http://de.mathworks.com/help/matlab/ref/fft.html
% FFT is calculated for binary sequences
function FFT                    % subfunction 3
load('S_bil_pups_binary_fs200.mat')
[~,m]=size(binarydata);  %dimensions of matrix n=number of rows, m=number of columns
gesamt_P1 = [];
gesamt_ioi =[];
for i= 1:m
k = find(binarydata(:,i)); % find outputs all values unequal zero -> finds all 1's, saves indices in k
X= binarydata(min(k):max(k),i); % X is defined as the binary sequence between the first 1 and the last 1 in the sequence
L = length(X); %L is the actual length of X from first onset to last onset, with deleted zeros at beginning and end

%get all iois
ioi=[];
for z=1:(length(k)-1)
    x= z+1;
    ioi(z,1)=k(x,1)-k(z,1);
end 


%% new Version
% https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/

X = 1/L*fftshift(fft(X,L)); %N-point complex DFT
df=FS/L;                    %frequency resolution
sampleIndex = -L/2:L/2-1;   %ordered index for FFT plot
f=sampleIndex*df;           %x-axis index converted to ordered frequencies

[pk1,lc1] = findpeaks(abs(X),'SortStr','descend','NPeaks',5); %finds five highest peaks
%highest peak in this is somewhere real close to 0 Hz -> Nyquist Frequency
P1peakP = X(lc1);           %Amplitude of five highest peaks
P1peakFreq = f(lc1);        %frequency of five highest peaks,
GOF = abs(P1peakP(2))/abs(P1peakP(1)); % goodness of fit (GOF) value, 
% we divide the amplitude of the second highest
% peak by the amplitude of the Nyquist Frequency which serves as an internal
% standard for what the maximum amplitude could actually be in the signal

%save all important information
bestf(i,1) = num2cell(abs(P1peakP(2)));         % FFT Amplitude
bestf(i,2) = num2cell(abs(P1peakFreq(2)));      % best fitting Pulse
bestf(i,3) = num2cell(L);                       % FFT Length
bestf(i,5) = num2cell(GOF);                     % Goodness of Fit
bestf(i,6)= {ioi};                              % IOI Sequenze for Model
end

%% old Version

% Y = fft(X);
% 
% %Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% % P1 based on P2 and the even-valued signal length L.
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% 
% 
% f = FS*(0:(L/2))/L; % 0-100 Hz, Original
% f=f'
% % OUTPUT: List of strongest responses with corresponding frequency
% % (descending order to easily get "best" pulse frequency
% Optimum=[];
% Optimum (:,1) = P1;
% Optimum (:,2) = f;
% % %Delet all values for f < 5
% %  indices = find(abs(Optimum(:,2)<5));
% %     Optimum(indices,:) = NaN; 
% %Delet all values for f < 1
%  indices = find(abs(Optimum(:,2)<1));
%     Optimum(indices,:) = NaN;
%  %find f for best(highest) P1
% [maxP1,ind] = max(Optimum(:,1));
% bestf(i,1) = num2cell(Optimum(ind,1));
% bestf(i,2) = num2cell(Optimum(ind,2));
% bestf(i,3) = num2cell(L); % amplitude should be correlated with L, the higher L the lower max amplitude
% bestf(i,5)= {ioi}


% save all amplitudes of all sequenzes to compare all against the maximal
% amplitudes found (histogram)
%gesamt_P1(end+1:end+size(P1))= P1;

gesamt_ioi(end+1:end+size(ioi,1),:) = ioi;


end 
listOfFileNames = listOfFileNames';
bestf(:,4)=listOfFileNames;             % original File Name

save('FFT_S_bil_pups_fs200_new.mat', ['bestf'])
% % to plot with stem 
% stem(f,abs(X));             %magnitudes vs frequencies
% xlabel('f (Hz)'); ylabel('|X(k)|');
end                             % subfunction 3



% %% 04.2 GOODNESS OF FIT
% % goodness of fit value is calculated: how? deviation of P from original
% % binary to "perfect" binary of that rhythm, same length, same element
% % number
% function gof_original           % subfunction 4
%     P_r_dev = {};
%    for n= 1: 500
%     X = []; %vector for 'perfect' binary sequence
%     R = 1000/cell2mat(bestf(n,2))/5; %value to produce sequence of perfect rhythm from original
%     L = cell2mat(bestf(n,3));
%     
%     X (1:L)= 0;
%     X (1)  = 1;
%     X (R:R:end)=1;
%         
%     k = find(X); % find outputs all values unequal zero -> finds all 1's, saves indices in k
%     X = X(min(k):max(k)); % X is defined as the binary sequence between the first 1 and the last 1 in the sequence
%     L = length(X);
%     
%     % jetzt FFT rechnen
%     Y = fft(X);
% 
% %Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% % P1 based on P2 and the even-valued signal length L.
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% 
% 
% f = FS*(0:(L/2))/L; % 0-100 Hz, Original
% f=f';
% % OUTPUT: List of strongest responses with corresponding frequency
% % (descending order to easily get "best" pulse frequency
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % change to "findpeaks" might be a lot easier!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Optimum=[];
% Optimum (:,1) = P1;
% Optimum (:,2) = f;
% % %Delet all values for f < 5
% %  indices = find(abs(Optimum(:,2)<5));
% %     Optimum(indices,:) = NaN; 
% %Delet all values for f < 1
%  indices = find(abs(Optimum(:,2)<1));
%     Optimum(indices,:) = NaN;
%  %find f for best(highest) P1
% [maxP1,ind] = max(Optimum(:,1));
% bestf_per(n,1) = num2cell(Optimum(ind,1));
% bestf_per(n,2) = num2cell(Optimum(ind,2));
%     % dann P_r_dev (P relativ deviation) berechnen
%     
%     P_r_dev{n,1} = abs(cell2mat(bestf_per(n,1))-cell2mat(bestf(n,1)))/cell2mat(bestf_per(n,1));
%     
%    end
%     
%    histogram(cell2mat(P_r_dev(:,1))) % show distribution of P_r_dev values close to 0 -> good, values close to 1 -> bad
%        
% end                             % subfcuntion 4    


%% 05.1 FFT MODEL
% random sequences are produced to test against original sequences
function random_sequence            %subfunction 5
%used data:
%bestf: from 'fft'
%gesamt_ioi: from 'fft'
binary_m=[]
binary_gesamt=[]

for seq_n=1:length(bestf)
    
for seq_n_i=1:50 %we want to produce 50 random sequences with the same element number as the original first sequence
    
    seqIOI=[];
    
for h=1:length(bestf{1,5}(:,1)) %bestf(~,5) contains the seqeunce of iois of the original seqeuence, we take the first sequence here
y=datasample(gesamt_ioi,1); %y is a random sample from gesamt_ioi, in bracktes: data to draw from(gesamt_ioi), number of datapoints to draw
seqIOI(h,1) = y; % is saved at corresponding step of sequence
end

%we need to convert the random iois back to a binary seqeunce
% just do the whole thing backwards: add all iois up
%example: ioi seqeunce: [18,13,17,15,13,16,15,20,15]
%then result should be: 1: 18
                        %2: 18+13
                        %3: 18+13+17
                        %4: ..
 for i=0:length(seqIOI)
     ind(i+1,1)= sum(seqIOI(1:i,1));
 end 
     
 % now we need a binary sequence 
 binary_m(1:max(ind), seq_n_i)=0; %corresponding column is set up in the right length, all zeros
 %indeces refered to in ind need to be '1'
 for i=1:length(ind)
 k= ind(i,1);
 binary_m(k+1,seq_n_i)=1;
 end
  [~,m]=size(binary_m);
end
%binary_gesamt(:,end+1:end+m) = binary_m; %soll gleich alle 500*50 random sequences in binary_gesamt speichern, anscheinend gibt es einen Dimension missmatch sobald eine sequence nicht die selbe maximale sequence länge hat wie vorher...vielleicht länger als alle mit 0 füllen vorher? oder einzeln speichern?

 save(['fft_model_seq', num2str(seq_n)], ['binary_m']);
end
end                             %subfunction 5
function fft_model              %subfunction 6
gesamt_P1 = [];
gesamt_ioi =[];
for seq = 1:500
    load(['fft_model_seq' , num2str(seq) '.mat'])
    [n,m]=size(binary_m);  %dimensions of matrix n=number of rows, m=number of columns
for i= 1:m
k = find(binary_m(:,i)); % find outputs all values unequal zero -> finds all 1's, saves indices in k
X= binary_m(min(k):max(k),i); % X is defined as the binary sequence between the first 1 and the last 1 in the sequence
L = length(X); %L is the actual length of X from first onset to last onset, with deleted zeros at beginning and end

%get all iois
ioi=[];
for z=1:(length(k)-1)
    x= z+1;
    ioi(z,1)=k(x,1)-k(z,1);
end 


Y = fft(X);

%Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);



f = FS*(0:(L/2))/L; % 0-100 Hz, Original
f=f';
% OUTPUT: List of strongest responses with corresponding frequency
% (descending order to easily get "best" pulse frequency
Optimum=[];
Optimum (:,1) = P1;
Optimum (:,2) = f;
% %Delet all values for f < 5
%  indices = find(abs(Optimum(:,2)<5));
%     Optimum(indices,:) = NaN; 
%Delet all values for f < 1
 indices = find(abs(Optimum(:,2)<1));
    Optimum(indices,:) = NaN;
 %find f for best(highest) P1
[maxP1,ind] = max(Optimum(:,1));
bestf_m(i,1) = num2cell(Optimum(ind,1));
bestf_m(i,2) = num2cell(Optimum(ind,2));
bestf_m(i,3) = num2cell(L); % amplitude should be correlated with L, the higher L the lower max amplitude
bestf_m(i,5)= {ioi};
% save all amplitudes of all sequenzes to compare all against the maximal
% amplitudes found (histogram)
gesamt_P1(end+1:end+size(P1))= P1;

gesamt_ioi(end+1:end+size(ioi,1),:) = ioi;
end 
save(['results_fft_model_seq', num2str(seq)], ['bestf_m']);
end


for a=1:500
load(['results_fft_model_seq', num2str(a)], ['bestf_m']);
mean1 = mean(cell2mat(bestf_m(1:50,1)));
gesamt_mean(a,1)= mean1;
end
save('mean_fft_p_pups','gesamt_mean')
end                                 %subfunction 6

% %% 05.2 FFT MODEL GOODNESS OF FIT
%     function gof_model          %subfunction 7
%         
%    P_r_dev_m = {};
%    
%    for m = 1:500
%        
%    load(['results_fft_model_seq', num2str(a)], ['bestf_m'])   
%    
%    for n= 1: 50
%     X = []; %vector for 'perfect' binary sequence
%     R = 1000/cell2mat(bestf_m(n,2))/5; %value to produce sequence of perfect rhythm from original
%     L = cell2mat(bestf_m(n,3));
%     
%     X (1:L)= 0;
%     X (1)  = 1;
%     X (R:R:end)=1;
%         
%     k = find(X); % find outputs all values unequal zero -> finds all 1's, saves indices in k
%     X = X(min(k):max(k)); % X is defined as the binary sequence between the first 1 and the last 1 in the sequence
%     L = length(X);
%     
%     % jetzt FFT rechnen
%     Y = fft(X);
% 
% %Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% % P1 based on P2 and the even-valued signal length L.
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% 
% 
% f = FS*(0:(L/2))/L; % 0-100 Hz, Original
% f=f';
% % OUTPUT: List of strongest responses with corresponding frequency
% % (descending order to easily get "best" pulse frequency
% Optimum=[];
% Optimum (:,1) = P1;
% Optimum (:,2) = f;
% % %Delet all values for f < 5
% %  indices = find(abs(Optimum(:,2)<5));
% %     Optimum(indices,:) = NaN; 
% %Delet all values for f < 1
%  indices = find(abs(Optimum(:,2)<1));
%     Optimum(indices,:) = NaN;
%  %find f for best(highest) P1
% [maxP1,ind] = max(Optimum(:,1));
% bestf_per_m(n,1) = num2cell(Optimum(ind,1));
% bestf_per_m(n,2) = num2cell(Optimum(ind,2));
%     % dann P_r_dev (P relativ deviation) berechnen
%     
%     P_r_dev_m{end+1,1} = abs(cell2mat(bestf_per_m(n,1))-cell2mat(bestf_m(n,1)))/cell2mat(bestf_per_m(n,1));
%         
%    end  
%    end  
%     end                         %subfunction 7


end 
