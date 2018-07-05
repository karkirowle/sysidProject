% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018

% Server progress
% Monitors progress of simulations

% Housekeeping
clc;
clear;
close all;

fileList = {'checkpoints/run_14_Jun_2018_07_11_15_50_200.mat'};

fileTimes = cell(1,3);

% Load progress data
load('serverprogress');


for f=1:size(fileList,2)
    % There are a huge amount of unimportant warnings due to simLib
    % objects, which we don't care about
    warning('off','all')
    load(fileList{f});
    
    % But let's be safe, so turn it back on
    warning('on','all');
    
    % Calculate remaining realisations - each SNR has same amount
    % of realisations
    totalRealisations = numRealisations * length(SNR);
    
    % How much is undone in current, minus and the remaining
    remainingRealisations = numRealisations* (length(SNR) - i) + ...,
        (numRealisations - r);
        
    % Percentage figures
    percentage =  100 - remainingRealisations./totalRealisations * 100;
    remaining = remainingRealisations./totalRealisations * 100;

    % Display filename and percentage
    disp(['Simulation Name: ', fileList{f}, newline, 'Progress: ', ...,
        num2str(percentage), '%']);
    
    
    prompt = 'Was this file downloaded recently? Y/N [N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    if strcmp(str,'Y')
        
        % Current cell size
        currentSize = size(fileTimes,1);
        
        % Check if file is in cell
        if any(strcmp(fileTimes,fileList{f}))
            id = strcmp(fileTimes,fileList{f});
            fileTimes{id,2} = clock;
            fileTimes{id,3} = remaining;
        else
            fileTimes{currentSize + 1,1} = fileList{f};
            fileTimes{currentSize + 1,2} = ...,
                [fileTimes{currentSize + 1,2},clock];
            fileTimes{currentSize + 1,3} = ...,
                [fileTimes{currentSize + 1,3},remaining];
            
            % At this branch it is possible to do ETA
            unwrapFileTimes = fileTimes{currentSize + 1,2};
            unwrapRemaining = fileTimes{currentSize + 1,2};
            
            timeDiff = etime(unwrapFileTimes(end-1), unwrapFileTimes(end));
            velocity = (unwrapRemaining(end-1)-unwrapRemaining)/timeDiff;
            ETA = remaining/velocity;
            
            
        end
        
    end
    
end

clearvars -except fileTimes
save('serverprogress')