function [data] = getEphSpaceTrackZip(files, StarlinkID, verbosity)
% getEphSpaceTrackZip gets the ephemeris data from the zip files provided
%                     by space-track.org
%
% -- Input --
% files            Nx1 cell of filenames corresponding to the space-track.org
%                  downloads. Each file should have a format such as:
%        
%                  SpaceX_Ephemeris_552_SpaceX_2022-12-12UTC21_21_02_1.zip
%        
%                  The most N can be is 3. The three files should only differ 
%                  by the last number (i.e.  ..._1.zip).
%
% StarlinkID       ID of starlink SV to grab
%
% verbosity        (optional) if set to anything will cause 
%                  command line output describing steps the function takes
%
% -- Output -- 
% A struct with the following fields:
%
% date_created       Datetime the ephemeris data was created in UTC timezone
%
% ephemeris_start    Datetime of the first ephemeris data point 
%
% ephemeris_stop     Datetime of the last ephemeris data point
%
% step_size          Length between ephemeris data points in seconds
%
% epoch_datetime     1x4321 vector of dates of specific epoch in UTC
%
% epoch_state        6x4321 vector of states of specific epoch. The state is
%                    defined as [x,y,z,vx,vy,vz]' in ECI.
%
% epoch_covariance   6x6x4321 vector of covariances of specific epoch state.

% Make sure file names are properly related
if (length(files)>3)
    error("Will only accept up to 3 files.")
end
filenames = cell(length(files),1);
directories = cell(length(files),1);
for ii = 1:length(files)
    [directories{ii},filenames{ii},ext] = fileparts(files{ii});
    if (string(ext) ~= ".zip")
        error("One or more of the files provided are not '.zip' files.")
    end
end
if (length(filenames) == 2)
    if (~strncmpi(filenames{1},filenames{2},length(files{1})-1))
        error("First two file names do not differ only by last value.")
    end
end
if (length(filenames) == 3)
    if (~strncmpi(filenames{1},filenames{2},length(filenames{1})-1))
        error("First two file names do not differ only by last value.")
    end
    if (~strncmpi(filenames{2},filenames{3},length(filenames{1})-1))
        error("The second and third file names do not differ only by last value.")
    end
    if (~strncmpi(filenames{1},filenames{3},length(filenames{1})-1))
        error("The first and third file names do not differ only by last value.")
    end
end
if ~exist('verbosity','var')
     % third parameter does not exist, so default it to something
     verbose = 0;
else 
    verbose = verbosity;
end 

pattern = sprintf("_STARLINK-%d_", StarlinkID);
if (verbose)
fprintf("Looking for STARLINK-%d...\n",StarlinkID);
end
fileWithId = {};
% Find zipped ephemeris
for ii = 1:length(files)
    if (~isempty(fileWithId))
        break
    end
    % Create a Java file of the ZIP filename.
    zipJavaFile  = java.io.File(files{ii});
    % Create a Java ZipFile and validate it.
    zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);
    % Extract the entries from the ZipFile.
    entries = zipFile.getEntries;
    % Loop through the entries and add to the file list.
    while entries.hasMoreElements
        txtFile = char(entries.nextElement);
        entry = zipFile.getEntry(txtFile); % empty if not found
        if (contains(string(txtFile),pattern))
            fileWithId = files{ii};
            % Create the Java File output object using the entry's name.
            outFile = strcat(string(directories{ii}),filesep,txtFile);
            if (verbose)
            fprintf("Found!\nExtracting to %s.\n",outFile);
            end
            inputstream = zipFile.getInputStream(entry);
            outJavaFile = java.io.File(fullfile(pwd,outFile));
            outStream = java.io.FileOutputStream(outJavaFile);
            copier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
            copier.copyStream(inputstream, outStream);
            outStream.close();
            zipFile.close
            if (verbose)
            fprintf("Success! Reading ephemeris file...\n");
            end
            break
        end
    end
    zipFile.close
end

if (isempty(fileWithId))
    warn = sprintf("STARLINK-%d not found in provided files.\n",StarlinkID);
    warning(warn);
    data = struct([]);
    return
end

% SpaceX uses Modified ITC format

pattern1 = sprintf("created:");
pattern2_1 = sprintf("ephemeris_start:");
pattern2_2 = sprintf("ephemeris_stop:");
pattern2_3 = sprintf("step_size:");
pattern3 = sprintf("ephemeris_source:");


fid = fopen(outFile,'r');
line = fgetl(fid);
iiEphemeris = 1;
HeaderEndedFlag = 0;
while ischar(line)
    if (contains(string(line),pattern1))
        data.date_created = datetime(erase(erase(string(line),pattern1),"UTC"),"TimeZone","UTC");
    end
    if (contains(string(line),pattern2_1))
        spl = split(line,pattern2_3);
        data.step_size = str2num(spl{2});
        spl = split(spl{1},pattern2_2);
        data.ephemeris_stop = datetime(erase(string(spl{2}),"UTC"),"TimeZone","UTC");
        spl = split(spl{1},pattern2_1);
        data.ephemeris_start = datetime(erase(string(spl{2}),"UTC"),"TimeZone","UTC");
        HeaderEndedFlag = 1;
    end
    if (HeaderEndedFlag)
        values = str2double(split(line,' '));
        % Valid data to save
        if (sum(isnan(values)) == 0)
            % Every ephemeris point has 1 line for
            % epoch date and time and state vector
            dt = string(values(1));
            dt = split(dt,'.');
            yyyyDOYhhmmss = char(dt(1));
            Y = str2num(yyyyDOYhhmmss(1:4));
            yyyyDOYhhmmss = yyyyDOYhhmmss(5:end);
            S = str2num(yyyyDOYhhmmss(end-1:end));
            yyyyDOYhhmmss = yyyyDOYhhmmss(1:end-2);
            M = str2num(yyyyDOYhhmmss(end-1:end));
            yyyyDOYhhmmss = yyyyDOYhhmmss(1:end-2);
            H = str2num(yyyyDOYhhmmss(end-1:end));
            DOY = str2num(yyyyDOYhhmmss(1:end-2));
            [~,Mo,D,~,~,~] = datevec(datenum(Y,1,DOY));
            if (length(dt) > 1)
                date = sprintf("%d-%d-%d %d:%d:%d.%d",Y,Mo,D,H,M,S,dt(2));
            else
                date = sprintf("%d-%d-%d %d:%d:%d",Y,Mo,D,H,M,S);
            end
            data.epoch_datetime(iiEphemeris) = datetime(date);
            data.epoch_state(:,iiEphemeris) = values(2:end);
            lowerDiagonal = [];
            for cov = 1:3
                line = fgetl(fid);
                if ~(ischar(line))
                    error("File is not the same length as expected or has been corrupted.")
                end
                values = str2double(split(line,' '));
                if (~sum(isnan(values)) == 0)
                    error("File is not the same length as expected or has been corrupted.")
                end
                lowerDiagonal = [lowerDiagonal; values];
            end
            a = triu(ones(6));
            a(a > 0) = lowerDiagonal;
            data.epoch_covariance(:,:,iiEphemeris) = (a + a')./(eye(6)+1);
            iiEphemeris = iiEphemeris + 1;
        end
    end
    line = fgetl(fid);
end
fclose(fid);
end