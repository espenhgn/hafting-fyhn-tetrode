% AxWaveForm2_1('inputFile.txt')
%
% Calculates the mean waveform for the 4 channels for each cell.
% From the mean waveform the peak to baseline amplititude and the peak lobe
% width are calcualted. The time between the peak and the trough after the
% peak is also calculated. The amplitudes are in micro volt.
%
% A variety of waveform parameters have been used to separate these two classes, including 
% trough-to-peak time (Bartho´ et al., 2004; Mitchell et al., 2007), 
% the ratio of trough-to-peak amplitude (Andermann et al., 2004; Hasenstaub et al., 2005), 
% and the duration of the longer, positive peak (Bruno and Simons, 2002; Atencio and Schreiner, 2008).
%
% If you set the parameter p.includeCellType to 1, you have to specify the 
% cell type for each cell. Each cell type will get 
% its own colour in the plot of peak to trough vs half width. Cell types
% are: Grid (Normal Grid cell), Conj (Conjunctive cell), 
% Hdir (Head direction cell), Bord (Border cell) and Undf (Undefined cell)
%
% Colour code:
% Grid  Red
% Hdir  Cyan
% Conj  Green
% Undf  Grey
% Bord  Magenta
%
% Cells with no trough after the peak is marked with an open circle.
%
%                   --- Input file example ---
%
% Session k:\haftingf\10586\16040402
% Tetrode 1
% Unit 1
% Type Grid
% Unit 2
% Type Conj
% ---
% Session k:\haftingf\10586\16040403
% Session k:\haftingf\10586\16040404
% Cut k:\haftingf\10586\16040403a_4.cut
% Tetrode 4
% Unit 1
% Type Hdir
% Unit 2
% Type Undf
%
%
%
%
%
% Version 1.1
% 14. May 2007
% 
% Version 1.2       Added plot of peak to trough vs half width for all the 
% 03. Nov. 2008     cells listed in the input file
%
% Version 1.3       Added colour code for the different cells.
% 08. Dec. 2008
%
% Version 1.4       Added columns for the highest ampitude channel in the
% 15. Dec. 2009     output file.
%                   Optimized the code for speed and memory usage.
%                   Added cubic spline interpolation of the mean wave form
%                   for each channel to get a finer resoultion of the plot
%                   and the calculation done on them.
%                   Added wave type to the output file.
%
% Version 1.5       Added option to run the program on input files that
% 11. Feb. 2010     don't contain cell type.
%                   Added the rat number to the output image file names.
%
% Version 1.6       Adjusted the plot.
% 12. Feb, 2010
%
% Version 2.0       Added convertion of the amplitude values from bits to
% 31. Jan. 2012     micro volt. Removed loading of position data that are
%                   not in use.
%
% Version 2.1       Torkel Hafting - added firing rate to output file, and
% 20. Feb 2013      merged figures to one file.

%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2007-2012.
function AxWaveForm2_1(inputFile)


%__________________________________________________________________________
%
%                       Program parameters
%__________________________________________________________________________


% Set if the input file will contain cell type or not
% 0 = No cell type in input file
% 1 = Cell type in the input file
p.includeCellType = 0;

% Format for the images
% format = 0 -> skip plotting procedure
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ai (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
% format = 7 -> fig (Matlab figure)
p.imageFormat = 4;

% Name for the folder where the images will be stored. In addition the name
% of the input file will be used in the folder name. Example: in121314.txt
% will give a folder name gridnessImages_in121314
p.imageFolder = 'WaveFormImages';


%__________________________________________________________________________





% Check if the output directory for the images is present, if not make it.
dirInfo = dir;
found = 0;
for kk=1:size(dirInfo,1)
    if dirInfo(kk).isdir
        if strcmp(dirInfo(kk).name,p.imageFolder)
            found = 1;
        end
    end
end
if found==0
    mkdir(p.imageFolder);
end


% Generate output filename
namelength=(length(inputFile)-4); 
corename=inputFile(3:namelength);
outputFile = strcat('out-',corename,'.xls');

% Open the input file for reading
fid = fopen(inputFile,'r');
if fid == -1 % Couldn't open the file
    disp('ERROR:');
    disp('      Could not open the input file!')
    disp('      Please make sure that the file name and path are correct.');
    return;
end
fid2 = fopen(outputFile,'w');
% Write header to the output file
fprintf(fid2,'%s\t','Session(s)');
fprintf(fid2,'%s\t','Tetrode');
fprintf(fid2,'%s\t','Unit');
fprintf(fid2,'%s\t','# spikes');
fprintf(fid2,'%s\t','Rate (Hz)');

fprintf(fid2,'%s\t','Width Highest Channel');
fprintf(fid2,'%s\t','Peak2Baseline Highest Channel');
fprintf(fid2,'%s\t','Peak2TroughTime Highest Channel');
fprintf(fid2,'%s\t','Peak2TroughAmplitude Highest Channel');
fprintf(fid2,'%s\t','Wave Type Highest Channel (0 = Normal, 1 = No trough after peak)');

fprintf(fid2,'%s\t','Width Ch1');
fprintf(fid2,'%s\t','Width Ch2');
fprintf(fid2,'%s\t','Width Ch3');
fprintf(fid2,'%s\t','Width Ch4');
fprintf(fid2,'%s\t','Peak2Baseline Ch1');
fprintf(fid2,'%s\t','Peak2Baseline Ch2');
fprintf(fid2,'%s\t','Peak2Baseline Ch3');
fprintf(fid2,'%s\t','Peak2Baseline Ch4');
fprintf(fid2,'%s\t','Peak2TroughTime Ch1');
fprintf(fid2,'%s\t','Peak2TroughTime Ch2');
fprintf(fid2,'%s\t','Peak2TroughTime Ch3');
fprintf(fid2,'%s\t','Peak2TroughTime Ch4');
fprintf(fid2,'%s\t','Peak2TroughAmplitude Ch1');
fprintf(fid2,'%s\t','Peak2TroughAmplitude Ch2');
fprintf(fid2,'%s\t','Peak2TroughAmplitude Ch3');
fprintf(fid2,'%s\t','Peak2TroughAmplitude Ch4');

fprintf(fid2,'\n');

% Set the maximum number of cells possible to have in the input file
maxNumCells = 100;

% Allocate memory for the arrays
peak2TroughTime = zeros(maxNumCells,1);
peak2baseline = zeros(maxNumCells,1);
peakWidth = zeros(maxNumCells,1);
cellCounter = 0;

% 0 = Undefined
% 1 = Grid
% 2 = Head Direction
% 3 = Conjunctive
cellType = zeros(maxNumCells,1);
% 0 = Normal
% 1 = No real trough after the peak
waveType = zeros(maxNumCells,1);

% Read information from the input file
while ~feof(fid)
    % Flag indicating if we are combining sessions into one big session. 0 =
    % single session, 1 = combined sessions with joint cut file, 2 =
    % combined session with separate cut files.
    combined = 0;
    % May have up to 10 combined session and cut files
    Xsessions = cell(10,1);
    Xcut = cell(10,1);
    Xcounter = 0;
    channelGain = cell(10,1);
    
    % Read first line of input file
    str = fgetl(fid);
    if length(str)<7
        disp('ERROR:');
        disp('  Expected string Session in input file, but was not found!');
        disp('FAILED!');
        return
    end
    if ~strcmpi(str(1:7),'Session')
        disp('ERROR');
        disp('  Expected the string Session in input file, but was not found!');
        disp('FAILED!');
        return
    else
        session = str(9:end);
        % Read next line
        str = fgetl(fid);
    end
    
    if strcmpi(str(1:3),'Cut')
        combined = 2;
        disp('Combined sessions with separate cut files detected')
        % Name for first cut file
        cutFile = str(5:end);
        % Read next line
        str = fgetl(fid);
        % Read session and cut info as long as there are more in the input
        % file
        getSess = 1;
        while 1
            if getSess % Session or room info expected next
                if strcmpi(str(1:7),'Session')
                    Xcounter = Xcounter + 1;
                    Xsessions{Xcounter} = str(9:end);
                    getSess = 0;
                    str = fgetl(fid);
                else
                    % No more session
                    break
                end
            else % Cut info expected next
                if strcmpi(str(1:3),'cut')
                    Xcut{Xcounter} = str(5:end);
                    str = fgetl(fid);
                    getSess = 1;
                else
                    msgbox('Missing cut file!','Error');
                    return;
                end
            end
        end
    end
    
    while strcmpi(str(1:7),'Session')
        % Sessions will be combined
        combined = 1;
        Xcounter = Xcounter + 1;
        Xsessions(Xcounter) = {str(9:end)};
        str = fgetl(fid);
    end
    
    if combined==1
        disp('Combined sessions with joint cut file detected');
        % Get shareed cut file for the combined sessions
        if strcmpi(str(1:3),'cut')
            cutFile = str(5:end);
            str = fgetl(fid);
        else
            msgbox('Missing cut file!','Error');
            return;
        end
    end
    
    
    fprintf('%s%s\n','Read data for session ',session);
    
    
    setFileName = sprintf('%s.set',session);
    [channelGain{1}, adcFullscale, duration] = getChannelGain(setFileName);
    
    if combined==1 || combined==2
         for ii = 1:Xcounter
             fprintf('%s%s\n','Read data for session ',Xsessions{ii});

             setFileName = sprintf('%s.set',Xsessions{ii});
             channelGain{ii+1} = getChannelGain(setFileName);
         end
    end
    
    % Continue to read data from the input file
    while ~feof(fid)
        if strcmp(str(1:3),'---')
            % New session will follow, start over.
            break;
        end
        if length(str) > 7
            if strcmp(str(1:7),'Tetrode') || strcmp(str(1:7),'tetrode') || strcmp(str(1:7),'TETRODE')
                tetrode = sscanf(str,'%*s %u');
                disp('Read spike data');
                
                if combined == 0
                    cutFile = sprintf('%s_%u.cut',session,tetrode);
                    cut = getcut(cutFile);
                elseif combined == 1
                    cut = getcut(cutFile);
                else
                    % Read the separate cut files, and join them together
                    cut = getcut(cutFile);
                    for ii = 1:Xcounter
                        cutTemp = getcut(Xcut{ii});
                        cut = [cut; cutTemp];
                    end
                end
                
                
                % Get spike data
                datafile = sprintf('%s.%u',session,tetrode);
                [~,ch1,ch2,ch3,ch4] = getspikes(datafile);
                
                % Channel gains for the 4 channels on this tetrode
                channels = (tetrode-1)*4+1:(tetrode-1)*4+4;
                gain = channelGain{1}(channels);
                
                % Convert the waveform amplitude from bit to micro volt
                ch1 = bits2Voltage(ch1, gain(1), adcFullscale, 1);
                ch2 = bits2Voltage(ch2, gain(2), adcFullscale, 1);
                ch3 = bits2Voltage(ch3, gain(3), adcFullscale, 1);
                ch4 = bits2Voltage(ch4, gain(4), adcFullscale, 1);
                
                 if combined==1 || combined==2
                    for ii = 1:Xcounter
                        datafile = sprintf('%s.%u',Xsessions{ii},tetrode);
                        [~,c1,c2,c3,c4] = getspikes(datafile);
                        
                        % Channel gains for the 4 channels on this tetrode
                        gain = channelGain{1+ii}(channels);
                        
                        % Convert the waveform amplitude from bit to micro volt
                        c1 = bits2Voltage(c1, gain(1), adcFullscale, 1);
                        c2 = bits2Voltage(c2, gain(2), adcFullscale, 1);
                        c3 = bits2Voltage(c3, gain(3), adcFullscale, 1);
                        c4 = bits2Voltage(c4, gain(4), adcFullscale, 1);

                        % Combine the wave form array
                        N1 = size(ch1,1);
                        N2 = size(c1,1);
                        combC1 = zeros(N1+N2,50);
                        combC2 = zeros(N1+N2,50);
                        combC3 = zeros(N1+N2,50);
                        combC4 = zeros(N1+N2,50);
                        combC1(1:N1,:) = ch1(:,:);
                        combC1(N1+1:end,:) = c1(:,:);
                        combC2(1:N1,:) = ch2(:,:);
                        combC2(N1+1:end,:) = c2(:,:);
                        combC3(1:N1,:) = ch3(:,:);
                        combC3(N1+1:end,:) = c3(:,:);
                        combC4(1:N1,:) = ch4(:,:);
                        combC4(N1+1:end,:) = c4(:,:);
                        ch1 = combC1;
                        ch2 = combC2;
                        ch3 = combC3;
                        ch4 = combC4;
                    end
                end

                % Read next line
                str = fgetl(fid);

                while length(str) > 4 && (strcmp(str(1:4),'Unit')||strcmp(str(1:4),'unit')||strcmp(str(1:4),'UNIT'))
                    % Unit number
                    unit = sscanf(str,'%*s %u');
                    
                    % Read next line
                    str = fgetl(fid);
                    
                    if p.includeCellType == 1
                        if length(str)>4 && strcmpi(str(1:4),'Type')

                            % Determine cell type
                            if length(str) < 9
                                disp('Error: The cell type string in the input file is in wrong format');
                            end

                            typeStr = str(6:9);
                            if strcmpi(typeStr,'undf')
                                cellType(cellCounter) = 0;
                            elseif strcmpi(typeStr,'grid')
                                cellType(cellCounter) = 1;
                            elseif strcmpi(typeStr,'hdir')
                                cellType(cellCounter) = 2;
                            elseif strcmpi(typeStr,'conj')
                                cellType(cellCounter) = 3;
                            elseif strcmpi(typeStr, 'bord')
                                cellType(cellCounter) = 4;
                            else
                                disp('Error: Cell type have to be undf, grid, hdir or conj')
                                return
                            end

                            % Read next line
                            str = fgetl(fid);
                        else
                            disp('ERROR: You have to set the cell type after the unit number');
                            return
                        end
                    end
                    
           % Call the main function on this cell
                    [P2TT,P2B,PW,WT] = mainFunc(ch1,ch2,ch3,ch4,cut,session,tetrode,unit,fid2,p,duration);
                    
                    if ~isempty(P2TT)
                        cellCounter = cellCounter + 1;
                    
                        peak2TroughTime(cellCounter) = P2TT;
                        peak2baseline(cellCounter) = P2B;
                        peakWidth(cellCounter) = PW;
                        waveType(cellCounter) = WT;
                    end
                end
            end
        end
    end
end

numCells = cellCounter;

if p.imageFormat > 0

if p.includeCellType == 1
    figure(5)
    for cc = 1:numCells
        switch cellType(cc)
            case 0
                if waveType(cc) == 0
                    plot(cc,peak2TroughTime(cc),'.','color',[0.5, 0.5, 0.5],'MarkerSize',20)
                    hold on
                else
                    plot(cc,peak2TroughTime(cc),'o','color',[0.5, 0.5, 0.5])
                    hold on
                    plot(cc,peak2TroughTime(cc),'.','color',[0.5, 0.5, 0.5])
                end
            case 1
                if waveType(cc) == 0
                    plot(cc,peak2TroughTime(cc),'.r','MarkerSize',16)
                    hold on
                else
                    plot(cc,peak2TroughTime(cc),'or')
                    hold on
                    plot(cc,peak2TroughTime(cc),'.r')
                end
            case 2
                if waveType(cc) == 0
                    plot(cc,peak2TroughTime(cc),'.c','MarkerSize',16)
                    hold on
                else
                    plot(cc,peak2TroughTime(cc),'oc')
                    hold on
                    plot(cc,peak2TroughTime(cc),'.c')
                end
            case 3
                if waveType(cc) == 0
                    plot(cc,peak2TroughTime(cc),'.g','MarkerSize',16)
                    hold on
                else
                    plot(cc,peak2TroughTime(cc),'og')
                    hold on
                    plot(cc,peak2TroughTime(cc),'.g')
                end
            case 4
                if waveType(cc) == 0
                    plot(cc,peak2TroughTime(cc),'.m','MarkerSize',16)
                    hold on
                else
                    plot(cc,peak2TroughTime(cc),'om')
                    hold on
                    plot(cc,peak2TroughTime(cc),'.m')
                end
        end

    end
    hold off
    xlabel('Cell number')
    ylabel('Peak to trough latency');
    axis([0,numCells+0.5, 0, max(peak2TroughTime)+0.2*max(peak2TroughTime)])
    sInd = strfind(inputFile,'.');
    figFile = sprintf('%s%s%u%s%u%s',inputFile(1:sInd(end)-1),'_P2T_Latency');
    imageStore(figure(5),p.imageFormat,figFile,300)


    figure(6)
    for cc = 1:numCells
        switch cellType(cc)
            case 0
                if waveType(cc) == 0
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.','color',[0.7, 0.7, 0.7],'MarkerSize',16)
                    hold on
                else
                    plot(peak2baseline(cc), peak2TroughTime(cc),'o','color',[0.7, 0.7, 0.7])
                    hold on
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.','color',[0.7, 0.7, 0.7])
                end
            case 1
                if waveType(cc) == 0
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.r','MarkerSize',16)
                    hold on
                else
                    plot(peak2baseline(cc), peak2TroughTime(cc),'or')
                    hold on
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.r')
                end
            case 2
                if waveType(cc) == 0
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.c','MarkerSize',16)
                    hold on
                else
                    plot(peak2baseline(cc), peak2TroughTime(cc),'oc')
                    hold on
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.c')
                end
            case 3
                if waveType(cc) == 0
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.g','MarkerSize',16)
                    hold on
                else
                    plot(peak2baseline(cc), peak2TroughTime(cc),'og')
                    hold on
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.g')
                end
            case 4
                if waveType(cc) == 0
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.m','MarkerSize',16)
                    hold on
                else
                    plot(peak2baseline(cc), peak2TroughTime(cc),'om')
                    hold on
                    plot(peak2baseline(cc), peak2TroughTime(cc),'.m')
                end
        end

    end
    hold off
    xlabel('Peak to baseline amplitude')
    ylabel('Peak to trough latency');
    axis([0,max(peak2baseline)+0.2*max(peak2baseline), 0, max(peak2TroughTime)+0.2*max(peak2TroughTime)])
    figFile = sprintf('%s%s%u%s%u%s',inputFile(1:sInd(end)-1),'_P2T_Latency_vs_P2B');
    imageStore(figure(6),p.imageFormat,figFile,300)

    figure(7)
    for cc = 1:numCells
        switch cellType(cc)
            case 0
                if waveType(cc) == 0
                    plot(peakWidth(cc), peak2TroughTime(cc),'.','color',[0.7, 0.7, 0.7],'MarkerSize',16)
                    hold on
                else
                    plot(peakWidth(cc), peak2TroughTime(cc),'o','color',[0.7, 0.7, 0.7])
                    hold on
                    plot(peakWidth(cc), peak2TroughTime(cc),'.','color',[0.7, 0.7, 0.7])
                end
            case 1
                if waveType(cc) == 0
                    plot(peakWidth(cc), peak2TroughTime(cc),'.r','MarkerSize',16)
                    hold on
                else
                    plot(peakWidth(cc), peak2TroughTime(cc),'or')
                    hold on
                    plot(peakWidth(cc), peak2TroughTime(cc),'.r')
                end
            case 2
                if waveType(cc) == 0
                    plot(peakWidth(cc), peak2TroughTime(cc),'.c','MarkerSize',16)
                    hold on
                else
                    plot(peakWidth(cc), peak2TroughTime(cc),'oc')
                    hold on
                    plot(peakWidth(cc), peak2TroughTime(cc),'.c')
                end
            case 3
                if waveType(cc) == 0
                    plot(peakWidth(cc), peak2TroughTime(cc),'.g','MarkerSize',16)
                    hold on
                else
                    plot(peakWidth(cc), peak2TroughTime(cc),'og')
                    hold on
                    plot(peakWidth(cc), peak2TroughTime(cc),'.g')
                end
            case 4
                if waveType(cc) == 0
                    plot(peakWidth(cc), peak2TroughTime(cc),'.m','MarkerSize',16)
                    hold on
                else
                    plot(peakWidth(cc), peak2TroughTime(cc),'om')
                    hold on
                    plot(peakWidth(cc), peak2TroughTime(cc),'.m')
                end
        end

    end

    hold off
    xlabel('Peak width')
    ylabel('Peak to trough latency');
    axis([0,max(peakWidth)+0.2*max(peakWidth), 0, max(peak2TroughTime)+0.2*max(peak2TroughTime)])
    figFile = sprintf('%s%s%u%s%u%s',inputFile(1:sInd(end)-1),'_P2T_Latency_vs_PeakWidth');
    imageStore(figure(7),p.imageFormat,figFile,300)
end

    elseif p.imageFormat==0;
    disp('Parameter set to skip making figures - data saved in xls file.');
end %plot procedure

% Close the files
fclose('all');
% close('all')
disp('Finished analysing');









function [peak2TroughTime,peak2baseline,peakWidth,waveType] = mainFunc(ch1,ch2,ch3,ch4,cut,session,tetrode,unit,fid2,p,duration)

fprintf('%s%u%s%u\n','Analysing cell ',unit,' on tetrode ',tetrode)


% Indexes to the spikes for the current cell
spikeInd = find(cut==unit);

%torkel Feb 2013
numberofspikes=length(spikeInd);
firingrate=numberofspikes/duration;

if isempty(spikeInd)
    disp('No spikes, skipping cell')
    peak2TroughTime = [];
    peak2baseline = [];
    peakWidth = [];
    waveType = [];
    return;
end

ind = strfind(session,'\');
ind2 = strfind(session,'/');
if isempty(ind) && isempty(ind2)
    saveName = sprintf('%s%s%u%s%u',session,'_T',tetrode,'C',unit);
else
    if ~isempty(ind)
        sess = session(ind(end)+1:end);
        ratName = session(ind(end-1)+1:ind(end)-1);
        saveName = sprintf('%s%s%s%s%s%s%u%s%u',p.imageFolder,'\',ratName,'_',sess,'_T',tetrode,'C',unit);
    end
    if ~isempty(ind2)
        sess = session(ind2(end)+1:end);
        ratName = session(ind2(end-1)+1:ind2(end)-1);
        saveName = sprintf('%s%s%s%s%s%s%u%s%u',p.imageFolder,'/',ratName,'_',sess,'_T',tetrode,'C',unit);
    end
end


% Select the waveforms that belong to the current cell
ch1 = ch1(spikeInd,:);
ch2 = ch2(spikeInd,:);
ch3 = ch3(spikeInd,:);
ch4 = ch4(spikeInd,:);

% Calculate the mean wave form for the 4 electrodes
ch1Mean = mean(ch1);
ch2Mean = mean(ch2);
ch3Mean = mean(ch3);
ch4Mean = mean(ch4);

% Interpolate the mean wave signal to get a better resolution
timeResolution = 0.002;
xx = 0:timeResolution:1-timeResolution;
x = 0:0.02:0.98;
ch1Mean = interp1(x,ch1Mean,xx,'spline');
ch2Mean = interp1(x,ch2Mean,xx,'spline');
ch3Mean = interp1(x,ch3Mean,xx,'spline');
ch4Mean = interp1(x,ch4Mean,xx,'spline');

% Calculate the base line for each mean wave
baseLine1 = mean(ch1Mean);
baseLine2 = mean(ch2Mean);
baseLine3 = mean(ch3Mean);
baseLine4 = mean(ch4Mean);

% Calculate the peak amplitude and peak position for each mean wave form
[peakCh1,peakIndCh1] = max(ch1Mean);
[peakCh2,peakIndCh2] = max(ch2Mean);
[peakCh3,peakIndCh3] = max(ch3Mean);
[peakCh4,peakIndCh4] = max(ch4Mean);

% Calculate the trough amplitude and the trough position for each mean wave
% form. This is the trough after the peak!
[troughCh1,troughIndCh1] = min(ch1Mean(peakIndCh1:end));
[troughCh2,troughIndCh2] = min(ch2Mean(peakIndCh2:end));
[troughCh3,troughIndCh3] = min(ch3Mean(peakIndCh3:end));
[troughCh4,troughIndCh4] = min(ch4Mean(peakIndCh4:end));


% Calculate the peak to trough
peak2TroughCh1 = peakCh1 - troughCh1;
peak2TroughCh2 = peakCh2 - troughCh2;
peak2TroughCh3 = peakCh3 - troughCh3;
peak2TroughCh4 = peakCh4 - troughCh4;

% Set the trough indexes according to the whole waveform
troughIndCh1 = troughIndCh1 + peakIndCh1 - 1;
troughIndCh2 = troughIndCh2 + peakIndCh2 - 1;
troughIndCh3 = troughIndCh3 + peakIndCh3 - 1;
troughIndCh4 = troughIndCh4 + peakIndCh4 - 1;

% Calculate the time difference between peak and trough
peak2TroughTimeCh1 = (troughIndCh1-peakIndCh1) * timeResolution;
peak2TroughTimeCh2 = (troughIndCh2-peakIndCh2) * timeResolution;
peak2TroughTimeCh3 = (troughIndCh3-peakIndCh3) * timeResolution;
peak2TroughTimeCh4 = (troughIndCh4-peakIndCh4) * timeResolution;

% Calculate the peak to baseline amplitude
peak2baselineCh1 = peakCh1 - baseLine1;
peak2baselineCh2 = peakCh2 - baseLine2;
peak2baselineCh3 = peakCh3 - baseLine3;
peak2baselineCh4 = peakCh4 - baseLine4;

% Calculate the half peak2baseline value
tresh1 = (peakCh1+baseLine1)/2;
tresh2 = (peakCh2+baseLine2)/2;
tresh3 = (peakCh3+baseLine3)/2;
tresh4 = (peakCh4+baseLine4)/2;




% Calculate the peak width
left = peakIndCh1;
right = peakIndCh1;
leftGo = 1;
rightGo = 1;
maxInd = length(ch1Mean);
while leftGo || rightGo
    % Check left
    if leftGo
        left = left - 1;
        if left==0
           left = 1;
           leftGo = 0;
        end
        if ch1Mean(left)<tresh1
            leftGo = 0;
            left = left + 1;
        end
    end
    % Check right
    if rightGo
        right = right + 1;
        % Bondery check
        if right>maxInd
            right = maxInd;
            rightGo = 0;
        end
        if ch1Mean(right)<tresh1
            rightGo = 0;
            right = right-1;
        end
    end
end
% Half peak width in seconds
peakWidth1 = (right-left+1)*timeResolution;
treshX1 = (left:right);

left = peakIndCh2;
right = peakIndCh2;
leftGo = 1;
rightGo = 1;
while leftGo || rightGo
    % Check left
    if leftGo
        left = left - 1;
        if left==0
           left = 1;
           leftGo = 0;
        end
        if ch2Mean(left)<tresh2
            leftGo = 0;
            left = left + 1;
        end
    end
    % Check right
    if rightGo
        right = right + 1;
        % Bondery check
        if right > maxInd
            right = maxInd;
            rightGo = 0;
        end
        if ch2Mean(right)<tresh2
            rightGo = 0;
            right = right-1;
        end
    end
end
% Half peak width in seconds
peakWidth2 = (right-left+1)*timeResolution;
treshX2 = (left:right);

left = peakIndCh3;
right = peakIndCh3;
leftGo = 1;
rightGo = 1;
while leftGo || rightGo
    % Check left
    if leftGo
        left = left - 1;
        if left==0
           left = 1;
           leftGo = 0;
        end
        if ch3Mean(left)<tresh3
            leftGo = 0;
            left = left + 1;
        end
    end
    % Check right
    if rightGo
        right = right + 1;
        % Bondery check
        if right > maxInd
            right = maxInd;
            rightGo = 0;
        end
        if ch3Mean(right)<tresh3
            rightGo = 0;
            right = right-1;
        end
    end
end
% Half peak width in seconds
peakWidth3 = (right-left+1)*timeResolution;
treshX3 = (left:right);

left = peakIndCh4;
right = peakIndCh4;
leftGo = 1;
rightGo = 1;
while leftGo || rightGo
    % Check left
    if leftGo
        left = left - 1;
        if left==0
           left = 1;
           leftGo = 0;
        end
        if ch4Mean(left)<tresh4
            leftGo = 0;
            left = left + 1;
        end
    end
    % Check right
    if rightGo
        right = right + 1;
        % Bondery check
        if right > maxInd
            right = maxInd;
            rightGo = 0;
        end
        if ch4Mean(right)<tresh4
            rightGo = 0;
            right = right-1;
        end
    end
end
% Half peak width in seconds
peakWidth4 = (right-left+1)*timeResolution;
treshX4 = (left:right);

% Set the values based on the larges peak to trough channel
[~,ind] = max([peak2TroughCh1,peak2TroughCh2,peak2TroughCh3,peak2TroughCh4]); 
ind = ind(1);
switch ind
    case 1
        peak2TroughTime = peak2TroughTimeCh1;
        peak2baseline = peak2baselineCh1;
        peakWidth = peakWidth1;
        % Set wave type
        if troughIndCh1 > length(ch1Mean) * 0.8
            waveType = 1;
        else
            waveType = 0;
        end
    case 2
        peak2TroughTime = peak2TroughTimeCh2;
        peak2baseline = peak2baselineCh2;
        peakWidth = peakWidth2;
        % Set wave type
        if troughIndCh2 > length(ch2Mean) * 0.8
            waveType = 1;
        else
            waveType = 0;
        end
    case 3
        peak2TroughTime = peak2TroughTimeCh3;
        peak2baseline = peak2baselineCh3;
        peakWidth = peakWidth3;
        % Set wave type
        if troughIndCh3 > length(ch3Mean) * 0.8
            waveType = 1;
        else
            waveType = 0;
        end
    case 4
        peak2TroughTime = peak2TroughTimeCh4;
        peak2baseline = peak2baselineCh4;
        peakWidth = peakWidth4;
        % Set wave type
        if troughIndCh4 > length(ch4Mean) * 0.8
            waveType = 1;
        else
            waveType = 0;
        end
end


% Find highest channel
tempArray = [peak2TroughCh1, peak2TroughCh2, peak2TroughCh3, peak2TroughCh4];
[~, ind] = max(tempArray);

% Write results to the output file
fprintf(fid2,'%s\t',session);
fprintf(fid2,'%u\t',tetrode);
fprintf(fid2,'%u\t',unit);
fprintf(fid2,'%u\t',numberofspikes);
fprintf(fid2,'%1.2f\t',firingrate);
switch ind(1)
    case 1
        fprintf(fid2,'%1.2f\t',peakWidth1);
        fprintf(fid2,'%3.2f\t',peak2baselineCh1);
        fprintf(fid2,'%1.3f\t',peak2TroughTimeCh1);
        fprintf(fid2,'%1.3f\t',peak2TroughCh1);
    case 2
        fprintf(fid2,'%1.2f\t',peakWidth2);
        fprintf(fid2,'%3.2f\t',peak2baselineCh2);
        fprintf(fid2,'%1.3f\t',peak2TroughTimeCh2);
        fprintf(fid2,'%1.3f\t',peak2TroughCh2);
    case 3
        fprintf(fid2,'%1.2f\t',peakWidth3);
        fprintf(fid2,'%3.2f\t',peak2baselineCh3);
        fprintf(fid2,'%1.3f\t',peak2TroughTimeCh3);
        fprintf(fid2,'%1.3f\t',peak2TroughCh3);
    case 4
        fprintf(fid2,'%1.2f\t',peakWidth4);
        fprintf(fid2,'%3.2f\t',peak2baselineCh4);
        fprintf(fid2,'%1.3f\t',peak2TroughTimeCh4);
        fprintf(fid2,'%1.3f\t',peak2TroughCh4);
end
fprintf(fid2,'%u\t',waveType);
fprintf(fid2,'%1.2f\t',peakWidth1);
fprintf(fid2,'%1.2f\t',peakWidth2);
fprintf(fid2,'%1.2f\t',peakWidth3);
fprintf(fid2,'%1.2f\t',peakWidth4);
fprintf(fid2,'%3.2f\t',peak2baselineCh1);
fprintf(fid2,'%3.2f\t',peak2baselineCh2);
fprintf(fid2,'%3.2f\t',peak2baselineCh3);
fprintf(fid2,'%3.2f\t',peak2baselineCh4);
fprintf(fid2,'%1.3f\t',peak2TroughTimeCh1);
fprintf(fid2,'%1.3f\t',peak2TroughTimeCh2);
fprintf(fid2,'%1.3f\t',peak2TroughTimeCh3);
fprintf(fid2,'%1.3f\t',peak2TroughTimeCh4);
fprintf(fid2,'%1.3f\t',peak2TroughCh1);
fprintf(fid2,'%1.3f\t',peak2TroughCh2);
fprintf(fid2,'%1.3f\t',peak2TroughCh3);
fprintf(fid2,'%1.3f\t',peak2TroughCh4);
fprintf(fid2,'\n');



% Session name
ind = strfind(session,'\');
sess = session(ind(end)+1:end);

if p.imageFormat >0; % Make plots of the mean wave form
    x = 0:timeResolution:1-timeResolution;
    figure
    subplot(2,2,1); plot(x,ch1Mean);
hold on
y = repmat(baseLine1,maxInd,1);
plot(x,y,'color',[0.5 0.5 0.5])
treshX1 = treshX1*timeResolution;
y = repmat(tresh1,length(treshX1));
plot(treshX1,y,'r')
hold off
xlabel('Wave duration [sec]')
ylabel('Amplitude [micro volt]');
titleString = sprintf('%s%u%s%u','Mean Waveform Channel 1 for cell T',tetrode,'C',unit);
    title(titleString)
    set(gcf,'Color',[1 1 1]);
%     figFile = sprintf('%s%s',saveName,'_Mean_Waveform_Ch1');
%     imageStore(figure(1),p.imageFormat,figFile,300)


% figure(2)
subplot(2,2,3);
plot(x,ch2Mean);
hold on
y = repmat(baseLine2,maxInd,1);
plot(x,y,'color',[0.5 0.5 0.5])
treshX2 = treshX2*timeResolution;
y = repmat(tresh2,length(treshX2));
plot(treshX2,y,'r')
hold off
xlabel('Wave duration [sec]')
ylabel('Amplitude [micro volt]');
titleString = sprintf('%s%u%s%u','Mean Waveform Channel 2 for cell T',tetrode,'C',unit);
title(titleString)
set(gcf,'Color',[1 1 1]);
% figFile = sprintf('%s%s',saveName,'_Mean_Waveform_Ch2');
% imageStore(figure(2),p.imageFormat,figFile,300)

% figure(3)
subplot(2,2,2);
plot(x,ch3Mean);
hold on
y = repmat(baseLine3,maxInd,1);
plot(x,y,'color',[0.5 0.5 0.5])
treshX3 = treshX3*timeResolution;
y = repmat(tresh3,length(treshX3));
plot(treshX3,y,'r')
hold off
xlabel('Wave duration [sec]')
ylabel('Amplitude [micro volt]');
titleString = sprintf('%s%u%s%u','Mean Waveform Channel 3 for cell T',tetrode,'C',unit);
title(titleString)
set(gcf,'Color',[1 1 1]);
% figFile = sprintf('%s%s',saveName,'_Mean_Waveform_Ch3');
% imageStore(figure(3),p.imageFormat,figFile,300)

% figure(4)
subplot(2,2,4);
plot(x,ch4Mean);
hold on
y = repmat(baseLine4,maxInd,1);
plot(x,y,'color',[0.5 0.5 0.5])
treshX4 = treshX4*timeResolution;
y = repmat(tresh4,length(treshX4));
plot(treshX4,y,'r')
hold off
xlabel('Wave duration [sec]')
ylabel('Amplitude [micro volt]');
titleString = sprintf('%s%u%s%u','Mean Waveform Channel 4 for cell T',tetrode,'C',unit);
title(titleString)
set(gcf,'Color',[1 1 1]);
% figFile = sprintf('%s%s',saveName,'_Mean_Waveform_Ch4');
figFile = sprintf('%s%s',saveName);
imageStore(gcf,p.imageFormat,figFile,300)
end





% Gets the adcFullScale value and the channel gain values from the setup
% file.
function [channelGain, adcFullscale, duration] = getChannelGain(setFileName)

channelGain = NaN(128,1);
adcFullscale = NaN;

% Open the setup file for reading
try
    fid = fopen(setFileName,'r');
catch me
    % Issue the Matlab error message
    throw(me);
end

% Make room for up to 10000 lines in the file.
setFile = cell(10000,1);
lineCount = 0;

% Read each line of the setfile and put it into the cell array
while ~feof(fid)
    lineCount = lineCount + 1;
    setFile{lineCount} = fgetl(fid);  
end

% Close the setup file
fclose(fid);

% Shorten the array to the length of the file
setFile = setFile(1:lineCount);

% torkel
triallength=setFile{5};
duration=str2num(triallength(9:end));

% Set the keyword to search for when finding the ADC fullscale value
keyword = 'ADC_fullscale_mv';
N = length(keyword);
for ii = 1:lineCount
    if length(setFile{ii}) > N
        if strcmpi(setFile{ii}(1:N), keyword)
            % ADC fullscale value found
            adcFullscale = str2double(setFile{ii}(N+2:end));
            break;
        end
    end
end


if adcFullscale ~= 3680 && adcFullscale ~= 1500
    disp('Error: Corrupt setup file. ADC fullscale value not found or recognized')
    return
end


keyword = 'gain_ch_';
N = length(keyword);

for ii = 1:lineCount
    if length(setFile{ii}) > N
        if strcmpi(setFile{ii}(1:N),keyword)
            % Found a channel gain line
            sInd = strfind(setFile{ii},' ');
            channelNumber = str2double(setFile{ii}(N+1:sInd-1));
            try
                channelGain(channelNumber+1) = str2double(setFile{ii}(sInd+1:end));
            catch me
                disp(me)
            end
        end
    end
end





% Transform the Axona signal from bits to micro volt
function amplituded = bits2Voltage(amplituded, gain, adcFullscale, bytesPerSample)

bits = 8 * bytesPerSample - 1;

ind = find(amplituded < 0);
amplituded(ind) = 1000 * (amplituded(ind) / (2^bits)) * (adcFullscale / gain);

ind = find(amplituded >= 0);
amplituded(ind) = 1000 * (amplituded(ind) / ((2^bits)-1)) * (adcFullscale / gain);






% Function for storing figures to file
% figHanle  Figure handle (Ex: figure(1))
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ai (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
% format = 7 -> fig (Matlab figure)
% figFile   Name (full path) for the file
% dpi       DPI setting for the image file
function imageStore(figHandle,format,figFile,dpi)

% Make the background of the figure white
set(figHandle,'color',[1 1 1]);
dpi = sprintf('%s%u','-r',dpi);

switch format
    case 1
        % Store the image as bmp (24 bit)
        figFile = strcat(figFile,'.bmp');
        print(figHandle, dpi, '-dbmp',figFile);
    case 2
        % Store image as png
        figFile = strcat(figFile,'.png');
        print(figHandle, dpi,'-dpng',figFile);
    case 3
        % Store image as eps (Vector format)
        figFile = strcat(figFile,'.eps');
        print(figHandle, dpi,'-depsc',figFile);
    case 4
        % Store image as jpg
        figFile = strcat(figFile,'.jpg');
        print(figHandle,dpi, '-djpeg',figFile);
    case 5
        % Store image as ai (Adobe Illustrator)
        figFile = strcat(figFile,'.ai');
        print(figHandle,dpi, '-dill',figFile);
    case 6
        % Store image as tiff (24 bit)
        figFile = strcat(figFile,'.tif');
        print(figHandle,dpi, '-dtiff',figFile);
    case 7
        % Store figure as Matlab figure
        figFile = strcat(figFile,'.fig');
        saveas(figHandle,figFile,'fig')
end


%__________________________________________________________________________
%
%                           Import rutines
%__________________________________________________________________________






function clust = getcut(cutfile)
fid = fopen(cutfile, 'rt');
clust = [];
while ~feof(fid)
    string = fgetl(fid);
    if ~isempty(string)
        if (string(1) == 'E') 
            break;
        end
    end
end
while ~feof(fid)
  string = fgetl(fid);
  if ~isempty(string)
     content = sscanf(string,'%u')';
     clust = [clust content];
  end
end
fclose(fid);
clust = clust';




function [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%
%   [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

[spikes,spikeparam] = importspikes(filename);
ts = [spikes.timestamp1]';
nspk = spikeparam.num_spikes;
spikelen = spikeparam.samples_per_spike;
ch1 = reshape([spikes.waveform1],spikelen,nspk)';
ch2 = reshape([spikes.waveform2],spikelen,nspk)';
ch3 = reshape([spikes.waveform3],spikelen,nspk)';
ch4 = reshape([spikes.waveform4],spikelen,nspk)';

function [spikes,spikeparam] = importspikes(filename)
%
%   [spikes,spikeparam] = importspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

fid = fopen(filename,'r','ieee-be');
if (fid < 0)
   error(sprintf('Could not open %s\n',filename)); 
end    

% read all bytes, look for 'data_start'
fseek(fid,0,-1);
sresult = 0;
[bytebuffer, bytecount] = fread(fid,inf,'uint8');
for ii = 10:length(bytebuffer)
    if strcmp( char(bytebuffer((ii-9):ii))', 'data_start' )
        sresult = 1;
        headeroffset = ii;
        break;
    end
end
if (~sresult)
    fclose(fid);
    error(sprintf('%s does not have a data_start marker', filename));
end

% count header lines
fseek(fid,0,-1);
headerlines = 0;
while(~feof(fid))
    txt = fgetl(fid);
    tmp = min(length(txt),10);
    if (length(txt))
        if (strcmp(txt(1:tmp),'data_start'))
            break;
        else
            headerlines = headerlines + 1;
        end
    else
        headerlines = headerlines + 1;
    end   
end    

% find timebase
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^timebase.*')))
        timebase = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 96 kHz');   
    timebase = 96000;    
end

% find duration
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^duration.*')))
        duration = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Duration not reported, defaulting to last time stamp');   
    duration = inf;    
end

% find number of spikes
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_spikes.*')))
        num_spikes = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of spikes not reported, using all that can be found');   
    num_spikes = inf;    
end

% find bytes per sample
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_sample.*')))
        bytes_per_sample = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per sample not reported, defaulting to 1');   
    bytes_per_sample = 1;    
end

% find bytes per timestamp
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_timestamp.*')))
        bytes_per_timestamp = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per timestamp not reported, defaulting to 4');   
    bytes_per_timestamp = 4;    
end

% find samples per spike
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^samples_per_spike.*')))
        samples_per_spike = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Samples per spike not reported, defaulting to 50');   
    samples_per_spike = 50;    
end

% check spike format
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^spike_format.*')))
        if (length(regexp(txt,'^spike_format t,ch1,t,ch2,t,ch3,t,ch4')))
            sresult = 1;
            break;
        else
           fclose(fid);
           error(sprintf('Unknown spike format, cannot read spikes from %s',filename));   
        end
    end
end    
if (~sresult)
    fclose(fid);
    error(sprintf('No spike format reported, cannot read spikes from %s.\nAre you sure this is a spike file?',filename));   
end

% close the file
fclose(fid);

% count the number of spikes in the file
spikelen = 4 * (bytes_per_sample * samples_per_spike + bytes_per_timestamp);
num_spikes_in_file = floor((bytecount - headeroffset)/spikelen);
if (isfinite(num_spikes))
    if (num_spikes_in_file > num_spikes)
        warning(sprintf('%d spikes reported in header, but %s seems to contain %d spikes.',num_spikes,filename,num_spikes_in_file));
    elseif (num_spikes_in_file < num_spikes)
        warning(sprintf('%d spikes reported in header, but %s can contain have %d spikes.',num_spikes,filename,num_spikes_in_file));
        num_spikes = num_spikes_in_file;    
    end
else
    num_spikes = num_spikes_in_file;
end
    
% allocate memory for return values

spikestruct = struct('timestamp1',0,'waveform1',zeros(samples_per_spike,1), ...
                     'timestamp2',0,'waveform2',zeros(samples_per_spike,1), ...
                     'timestamp3',0,'waveform3',zeros(samples_per_spike,1), ...
                     'timestamp4',0,'waveform4',zeros(samples_per_spike,1));

spikes = repmat(spikestruct,num_spikes,1);
                        
% out the spikes into the struct, one by one

big_endian_vector =  (256.^((bytes_per_timestamp-1):-1:0))';
little_endian_matrix = repmat(256.^(0:(bytes_per_sample-1))',1,samples_per_spike);

for ii = 1:num_spikes
   % sort the bytes for this spike
   spikeoffset = headeroffset + (ii-1)*spikelen;
   t1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w1_bytes( w1_bytes > 127 ) = w1_bytes( w1_bytes > 127 ) - 256;
   w1_bytes = reshape(w1_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w2_bytes( w2_bytes > 127 ) = w2_bytes( w2_bytes > 127 ) - 256;
   w2_bytes = reshape(w2_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w3_bytes( w3_bytes > 127 ) = w3_bytes( w3_bytes > 127 ) - 256;
   w3_bytes = reshape(w3_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w4_bytes( w4_bytes > 127 ) = w4_bytes( w4_bytes > 127 ) - 256;
   w4_bytes = reshape(w4_bytes,bytes_per_sample,samples_per_spike);
   % interpret the bytes for this spike
   spikes(ii).timestamp1 = sum(t1_bytes .* big_endian_vector) / timebase; % time stamps are big endian
   spikes(ii).timestamp2 = sum(t2_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp3 = sum(t3_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp4 = sum(t4_bytes .* big_endian_vector) / timebase;
   spikes(ii).waveform1 =  sum(w1_bytes .* little_endian_matrix, 1); % signals are little-endian
   spikes(ii).waveform2 =  sum(w2_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform3 =  sum(w3_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform4 =  sum(w4_bytes .* little_endian_matrix, 1);
end
if (~isfinite(duration))
    duration = ceil(spikes(end).timestamp1);
end
spikeparam = struct('timebase',timebase,'bytes_per_sample',bytes_per_sample,'samples_per_spike',samples_per_spike, ...
                    'bytes_per_timestamp',bytes_per_timestamp,'duration',duration,'num_spikes',num_spikes);


 

%__________________________________________________________________________
%
%           Function for fixing the position timestamps
%__________________________________________________________________________

function [didFix,fixedPost] = fixTimestamps(post)

% First time stamp in file
first = post(1);
% Number of timestamps
N = length(post);
uniqePost = unique(post);

if length(uniqePost)~=N
    didFix = 1;
    numZeros = 0;
    % Find the number of zeros at the end of the file
    while 1
        if post(end-numZeros)==0
            numZeros = numZeros + 1;
        else
            break;
        end
    end
    
    last = first + (N-1-numZeros) *0.02;
    fixedPost = first:0.02:last;
    fixedPost = fixedPost';
else
    didFix = 0;
    fixedPost = [];
end




