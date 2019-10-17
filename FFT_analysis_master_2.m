%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFT ANALYSIS
%author: Lara S. Burchardt, l.s.burchardt@gmx.de
%version: 0.2.0 
%28.01.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fft_analysis_master_2 %main function

SetPath                     %subfunction 1, labeled at beginning and end
binary                      %subfunction 2, labeled at beginning and end
FFT                         %subfunction 3, labeled at beginning and end

%% 01: GLOBALS
% definde global variables

global listOfFileNames

%% 02: DATA MANAGMENT
% load data
% Initialization steps:
function SetPath        %subfunction 1

workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;

% Step one: Define a starting folder.
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
%filePattern = sprintf('%s/**/*.xlsx', topLevelFolder); %for .xlsx data
filePattern = sprintf('%s/**/*.xls', topLevelFolder); %for .xls data
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

prompt = 'Assign a name for saving. Example: C_persp_IC for Isolation Calls of Carollia perspicillata. The name needs to be put in aposthrophes. ';
savename = input(prompt);
fprintf(['The name you chose for saving is ' savename]);
end                         % subfunction 1


%% 03: TRANSFORM TO BINARY
% timestamps are converted into a binary sequence

function binary             % subfunction 2

%Step three: loop through files, read .xls files (xlsread) and calculate binary
%code for syllable onsets

binarydata = [];
binarydata(:,length(listOfFileNames)) = 0;

for k= 1: length(listOfFileNames)
    

    matfilename = listOfFileNames{:,k};
   
    data = xlsread(matfilename); 
    vectorA = data (:,1); % we just need the syllable onsets 
    vectorA = round(vectorA * 1000 / 5); %we multiply by 1000 to get data 
    % in miliseconds, divide by 5/1 because that is the timeresolution i want
    % to have and round to integer (200 Hz sampling rate--> sample every 5 milliseconds, sampling rate 1000 Hz--> sample every milisecond)

    
    for l = 1:max(vectorA)
        for i = 1:length(vectorA)
            if l == vectorA(i, 1)
            binarydata(l,k) = 1 ; 
            %else binarydata (l,k) = 0;
            end 
        end 
    end
    
% binarydata = array2table(binarydata); %converts array of kind 'double' to datatabel
                                      % uses less memory, easier to
                                      % manipulate later
end

xlswrite(['Binary_' savename '_fs200.xlsx'], binarydata);
save(['Binary_' savename '_fs200.mat']);


end                             % subfunction 2



%% 04.1: FFT
% Step four: FFT is calculated for binary sequences
function FFT                    % subfunction 3
load(['Binary_' savename '_fs200.mat'], 'binarydata')
[~,m]=size(binarydata);  %dimensions of matrix n=number of rows, m=number of columns

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
FS = 200;
X = 1/L*fftshift(fft(X,L)); %N-point complex DFT
df=FS/L;                    %frequency resolution
sampleIndex = -L/2:L/2-1;   %ordered index for FFT plot
f=sampleIndex*df;           %x-axis index converted to ordered frequencies

[~,lc1] = findpeaks(abs(X),'SortStr','descend','NPeaks',10); %finds five highest peaks
%highest peak in this is somewhere real close to 0 Hz -> Nyquist Frequency
P1peakP = X(lc1);           %Amplitude of five highest peaks
P1peakFreq = f(lc1);        %frequency of five highest peaks

if P1peakFreq(1,1) ~=0                          %account for shift in zero-bin component
P1peakFreq(1,:) = P1peakFreq - P1peakFreq(1,1); % gets shifted back to 0 here
else 
P1peakFreq(1,:) = P1peakFreq;
end

%Step five: calculate goodness-of-fit values 
if length(P1peakP) > 1
    GOF     = abs(P1peakP(2))/abs(P1peakP(1)); 

% goodness of fit (GOF) value, 
% we divide the amplitude of the second highest
% peak by the amplitude of the Nyquist Frequency which serves as an internal
% standard for what the maximum amplitude could actually be in the signal
% version 2, normalize peak 2 for number of samples

    GOF_2   = abs(P1peakP(2))/(L*abs(P1peakP(1)));

else
    GOF     = NaN;
    GOF_2   = NaN;
end



%Step six: save all important information
if length(P1peakP) > 1
bestf(i,1) = num2cell(abs(P1peakP(2)));         % FFT Amplitude
bestf(i,2) = num2cell(abs(P1peakFreq(2)));      % best fitting Pulse
bestf(i,5) = num2cell(GOF);                     % Goodness of Fit
bestf(i,7) = num2cell(GOF_2);                   % Goodness of Fit, version 2
else
    bestf(i,1) = {NaN};                         
    bestf(i,2) = {NaN};
    bestf(i,5) = {NaN};
    bestf(i,7) = {NaN};
end 
    bestf(i,3) = num2cell(L);                    % FFT Length
    bestf(i,6)= {ioi};                           % IOI Sequenze for Model
end

%total_ioi(end+1:end+size(ioi,1),:) = ioi;

end                                             % subfunction 3

listOfFileNames = listOfFileNames';
bestf(:,4)=listOfFileNames;                     %original FileName


save(['FFT_' savename '_fs200.mat'], ['bestf']);    %save 'bestf' in .mat file
xlswrite(['FFT_' savename '_fs200.xlsx'], bestf);   %save 'bestf' in .xlsx file

% % to plot with stem 
% stem(f,abs(X));             %magnitudes vs frequencies
% xlabel('f (Hz)'); ylabel('|X(k)|');
end                             






