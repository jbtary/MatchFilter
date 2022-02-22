% Match filter code for the Sea of Marmara
% First processing step using cc on complete events
% Portion and stations for each parent event selected by picking

% Suggestions:
% Rotate data in P-wave frame of reference?

clear
% SAC routines
addpath('/home/jbtary/Documents/MATLAB/SAC')

fsmaster = 125; % Sampling frequency Hz

% Band-pass filtering
fNy = fsmaster/2;
fcl = 15; fch = 25; % Frequency bounds in Hz
dirgen = '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/';
ndir = dir([dirgen 'G*']); ndir2 = dir([dirgen 'K*']);
ndir = [ndir;ndir2];
nametmp = extractfield(ndir,'name'); clear ndir2
savepath = '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ChildEvents3/'; % Path to folder to save child events

[b,a]=butter(6,[fcl fch]/fNy,'bandpass'); % Butterworth filter of order 6

% Definition of parent events
infopar = dir([dirgen 'ParDir/Picked/*.mat']); % Inventory of parent events (from Cut_Pick_ParentEvents.m)
filele = 3600; % Normal length of continuous data files in seconds

% Detection params (see synth_mad.m to define them)
madthres = 10; % MAD multiplifier
ntime = 10*fsmaster; % Min time between consecutive peaks/events in samples
win = 10*fsmaster; % Time before and after cross-corr max to take with child events
winP = 125; winS = 190; % Windows around P and S-waves for cross-corr. in samples

% Get list of data files concatenated in one file for all stations
infoc = format_filelists(dirgen,datenum([2011,04,15,12,0,0]),datenum([2011,07,31,12,0,0]));

fid = fopen([savepath '/Process_MatchFilter.txt'],'w');
tracability = zeros(size(infoc,1)*size(infopar,1),3);
for ii = 1:size(infoc,1) % Loop on the continuous data files
    tic
    
    kk = 1;
    for jj = 1:size(infoc,2) % Get data file
        
        if cell2mat(infoc(ii,jj)) ~= 0
            % Check on KOERI or OBS data (3,4,2 or e,n,z)
            if strcmp(nametmp{jj}(1),'G') == 1
                fs = 125; % Sampling frequency Hz
                compo = [3 4 2];
                path1 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'X/' infoc{ii,jj}(1:23) num2str(compo(1)) '.sac'];
                filepath1 = dir(path1);
                path2 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'Y/' infoc{ii,jj}(1:23) num2str(compo(2)) '.sac'];
                filepath2 = dir(path2);
                path3 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'Z/' infoc{ii,jj}(1:23) num2str(compo(3)) '.sac'];
                filepath3 = dir(path3);
                
                % Check on missing files
                if isempty(filepath1) == 1; clear filepath1; filepath1.bytes = 0; end
                if isempty(filepath2) == 1; clear filepath2; filepath2.bytes = 0; end
                if isempty(filepath3) == 1; clear filepath3; filepath3.bytes = 0; end
                
                % Check on empty files
                check = [filepath1.bytes filepath2.bytes filepath3.bytes];
            else
                fs = 100; % Sampling frequency Hz
                compo = ['e';'n';'z'];
                path1 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'E/' infoc{ii,jj}(1:23) num2str(compo(1)) '.sac'];
                filepath1 = dir(path1);
                path2 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'N/' infoc{ii,jj}(1:23) num2str(compo(2)) '.sac'];
                filepath2 = dir(path2);
                path3 = [dirgen ndir(jj).name '/' ...
                    ndir(jj).name 'Z/' infoc{ii,jj}(1:23) num2str(compo(3)) '.sac'];
                filepath3 = dir(path3);
                
                if isempty(filepath1) == 1; clear filepath1; filepath1.bytes = 0; end
                if isempty(filepath2) == 1; clear filepath2; filepath2.bytes = 0; end
                if isempty(filepath3) == 1; clear filepath3; filepath3.bytes = 0; end
                
                check = [filepath1.bytes filepath2.bytes filepath3.bytes];
            end
            
            % Read files, if one is zero, replace it by zeros
            for rr = 1:3 % Barbaric way
                if check(rr) ~= 0
                    eval(['[data' num2str(rr) ',hd]=rdSac(path' num2str(rr) ');'])
                else
                    eval(['data' num2str(rr) '=zeros(filele*fs,1);'])
                end
            end
            clear filepath* path* check
            
            % Check on the length of data files
            if length(data1) ~= filele*fs || length(data2) ~= filele*fs || length(data3) ~= filele*fs
                disp(['File ' infoc{ii,jj}(1:23) ' wrong length.']);
                continue;
            end
            
            % Re-sample data to concatenate it (useful when different
            % data have different sampling frequency)
            if fs ~= fsmaster
                % /!\ be mindful of end-effects of resample
                data1 = resample(data1,fsmaster,fs);
                data2 = resample(data2,fsmaster,fs);
                data3 = resample(data3,fsmaster,fs);
            end
            
            datax(kk,:) = data1';
            datay(kk,:) = data2';
            dataz(kk,:) = data3';
            
            sta(kk) = jj; % Track of station indexes
            kk = kk + 1;
            clear data1 data2 data3 fs rr compo
        end
    end
    clear kk jj
    
    dataxc = datax; datayc = datay; datazc = dataz; % Saving raw data
    
    % Filter continuous data
    datax = detrend(datax'); % Detrend and transpose array (detrend and filtfilt work on columns)
    datax = filtfilt(b,a,datax); datax = datax';
    datay = detrend(datay');
    datay = filtfilt(b,a,datay); datay = datay';
    dataz = detrend(dataz');
    dataz = filtfilt(b,a,dataz); dataz = dataz';
    
    % Normalize by avg. envelope (from E. Caffagni MFA code)
    for yy = 1:size(datax,1)
        env1 = envsm(datax(yy,:),[1,2],30); % Cont. data
        env2 = envsm(datay(yy,:),[1,2],30);
        env3 = envsm(dataz(yy,:),[1,2],30);
        
        env = (env1+env2+env3)/3;
        datax(yy,:) = datax(yy,:)./env;
        datay(yy,:) = datay(yy,:)./env;
        dataz(yy,:) = dataz(yy,:)./env; clear env*
    end
    clear yy
    
    for pp = 1:size(infopar,1) % Loop on parent events
        
        % Data of parent event
        parent = load([dirgen 'ParDir/Picked/' infopar(pp).name]);
        lettre = char(parent.comp);
        
        % Get corresponding continuous data
        [~,ia,ib] = intersect(nametmp(sta),lettre);
        datmpx = datax(ia,:);
        datmpy = datay(ia,:);
        datmpz = dataz(ia,:);
        
        dataxs = parent.dataxs(ib,:); % Select matching data from parent event
        datays = parent.datays(ib,:); % These are cells still
        datazs = parent.datazs(ib,:);
                
        % T0 of picks
        t0 = datenum(parent.newhdr(9),parent.newhdr(10),parent.newhdr(11),...
            parent.newhdr(12),parent.newhdr(13),parent.newhdr(14));
        var = parent.ppicks; var2 = parent.stap; var3 = parent.spicks; var4 = parent.stas;
        
        stapc = []; stasc = [];
        for jj=1:size(var2,1); stapc{1,jj}=var2(jj,:); end;
        for jj=1:size(var4,1); stasc{1,jj}=var4(jj,:); end;
        
        kp = 1; ks = 1; tp = []; ts = [];
        stasvp = []; stasvs = [];
        tshiftc = NaN(size(dataxs,1),2); % Store the time shifts of P and S wave picks
        for yy = 1:size(dataxs,1) % Loop on number of stations
            
            dataxp = dataxs{yy,1}; datayp = datays{yy,1}; datazp = datazs{yy,1};
            % Filter parent data
            dataxp = detrend(dataxp'); dataxp = filtfilt(b,a,dataxp); dataxp = dataxp';
            datayp = detrend(datayp'); datayp = filtfilt(b,a,datayp); datayp = datayp';
            datazp = detrend(datazp'); datazp = filtfilt(b,a,datazp); datazp = datazp';
            
            env1 = envsm(dataxp,[1,2],30); % Parent event data
            env2 = envsm(datayp,[1,2],30);
            env3 = envsm(datazp,[1,2],30);
            env = (env1+env2+env3)/3;
            dataxp=dataxp./env; datayp=datayp./env; datazp=datazp./env; clear env*
            
            stacurrent = parent.comp{1,yy}; % Current station
            [~,~,ip] = intersect(stacurrent,stapc,'stable'); % Position of station in picks
            [~,~,is] = intersect(stacurrent,stasc,'stable'); % Position of station in picks
            % P-waves
            if isempty(ip) ~= 1
                varcurrent = var(ip);
                % Time to picks in cut datastream in samples:
                tb=round(((t0+(varcurrent/86400))-parent.datebeg(yy))*86400*fsmaster);
                te=tb+winP;
                % Selection of data segment after P and/or S wave
                if tb<1; tb=1; end; if te>length(dataxp); te=length(dataxp); end; 
                dataxpp = dataxp(:,tb:te);
                dataypp = datayp(:,tb:te);
                datazpp = datazp(:,tb:te);
                
                % Lag of data2 vs data1
                [xcx,~] = xcorr(datmpx(yy,:),dataxpp);
                [xcy,~] = xcorr(datmpy(yy,:),dataypp);
                [xcz,lag] = xcorr(datmpz(yy,:),datazpp);
                
                xcx = xcx(1,lag>=0); xcy = xcy(1,lag>=0); xcz = xcz(1,lag>=0);
                xsump(kp,:) = xcx + xcy + xcz; % Sum all comps
                xsump(kp,:) = xsump(kp,:)/std(xsump(kp,:));
                % Save absolute time picks for later
                tp(kp,1) = t0+(varcurrent/86400);
                stasvp{1,kp} = stacurrent; % Keep track of stations
                
                kp = kp+1; % Increment for next station
                
                tshiftc(yy,1) = t0+(varcurrent/86400);
                clear ip varcurrent tb te dataxpp dataypp datazpp xcx xcy
                clear xcz lag
            end
            
            % S-waves
            if isempty(is) ~= 1
                varcurrent = var3(is);
                % Time to picks in cut datastream in samples:
                tb=round(((t0+(varcurrent/86400))-parent.datebeg(yy))*86400*fsmaster);
                te=tb+winS;
                % Selection of data segment after P and/or S wave
                if tb<1; tb=1; end; if te>length(dataxp); te=length(dataxp); end; 
                dataxpp = dataxp(:,tb:te);
                dataypp = datayp(:,tb:te);
                datazpp = datazp(:,tb:te);
                
                % Lag of data2 vs data1
                [xcx,~] = xcorr(datmpx(yy,:),dataxpp);
                [xcy,~] = xcorr(datmpy(yy,:),dataypp);
                [xcz,lag] = xcorr(datmpz(yy,:),datazpp);
                
                xcx = xcx(1,lag>=0); xcy = xcy(1,lag>=0); xcz = xcz(1,lag>=0);
                xsums(ks,:) = xcx + xcy + xcz; % Sum all comps
                xsums(ks,:) = xsump(ks,:)/std(xsump(ks,:));
                % Save absolute time picks for later
                ts(ks,1) = t0+(varcurrent/86400);
                stasvs{1,ks} = stacurrent; % Keep track of stations
                
                ks = ks+1; % Increment for next station
                
                tshiftc(yy,2) = t0+(varcurrent/86400);
                clear is varcurrent tb te dataxpp dataypp datazpp xcx xcy
                clear xcz lag
            end
            
            clear dataxp datayp datazp stacurrent ip is
        end
        
        % Concat absolute times associated with cross-corr datastreams
        tcc = [tp;ts]; % Both P and S CC streams at the same time
        stasv = [stasvp stasvs];
        tshiftc = round((tshiftc - min(tcc))*86400*fsmaster); % To keep P and S wave measures and station/data structure
        
        % Shift time series before summing (to correct for moveout of event)
        tshift = round((tcc - min(tcc))*86400*fsmaster);
        xsumc = [xsump;xsums];
        xcnew = zeros(size(xsumc,1),size(xsumc,2)+max(tshift));
        
        for yy = 1:length(tshift)
            xcnew(yy,max(tshift)-tshift(yy)+1:max(tshift)-tshift(yy)+size(xsumc,2)) = xsumc(yy,:);
        end
        
        xsum = sum(xcnew,1); % Sum all stations
        xmad = mad(xsum,1); % Median absolute deviation of xsum
        xsum = xsum(:,max(tshift)+1:end);
        
        % Find where we have values > thres separated by ntime
        [peaks,locs] = findpeaks(xsum,'MINPEAKHEIGHT',...
            madthres*xmad,'MINPEAKDISTANCE',ntime);
        
        hdrpar = parent; % Removing unnecessary fields to save disk space
        hdrpar = rmfield(hdrpar,'dataxp'); hdrpar = rmfield(hdrpar,'datayp'); 
        hdrpar = rmfield(hdrpar,'datazp');
        hdrpar = rmfield(hdrpar,'dataxs'); hdrpar = rmfield(hdrpar,'datays');
        hdrpar = rmfield(hdrpar,'datazs');
        hdrpar = rmfield(hdrpar,'datebegp'); hdrpar = rmfield(hdrpar,'compp');
        
        % Save events
        dirname = infopar(pp).name;
        if isempty(locs) ~= 1
            % Create directory for current parent event if needed
            savepathtmp = [savepath dirname(1:end-4)];
            if exist(savepathtmp,'dir') ~= 7
                mkdir(savepathtmp);
            end
            
            % Reseting to raw data for saving
            datmpx = dataxc(ia,:); datmpy = datayc(ia,:); datmpz = datazc(ia,:);
            
            % Correct from moveout of the different stations, get only
            % one value per station (get preferably S-wave measures
            % because they usually have higher amplitudes)
            tshiftc=max(tshiftc,[],2);
            
            % Sometimes there is not any data for a station and
            % tshiftc=NaN -> should actually be updated so that picks
            % and data correspond exactly. Only corrected for now.
            hdr.comp = nametmp(sta); % Comps saved for all comps (ex: from dataxc)
            comp = hdr.comp(ia); % Comps corresponding to data in cells like dataxs
            [stab,~,~] = unique(stasv,'stable'); % Stations in tshiftc (non-sorted)
            % Get data for stations in both tshiftc(stab) and comp (dataxs)
            [compsv,ua,ub] = intersect(comp,stab,'stable');
            % Select data
            datmpx = datmpx(ua,:); datmpy = datmpy(ua,:); datmpz = datmpz(ua,:);
            % Make sure that compsv and tshiftc are in the same order as
            % data for cutting and saving
            [~,va,vb] = intersect(stab,compsv,'stable');
            stab = stab(vb);
            tshiftc = tshiftc(vb);
            clear ua ub va vb
            
            if length(stab) ~= length(tshiftc)
                disp('Call to the operator, stab and tshiftc dont have the same length. In pause.')
                pause;
            end
            
            if size(datmpx,1) ~= length(tshiftc)
                disp('Call to the operator, datmpx and tshiftc dont have the same sta #. In pause.')
                pause;
            end
                        
            % Get rid of data without any detection
            ind = find(isnan(tshiftc)==1);
            if isempty(ind)~=1
                tshiftc(ind) = [];
                compsv(ind) = [];
                datmpx(ind,:) = []; 
                datmpy(ind,:) = []; 
                datmpz(ind,:) = [];
            end
            
            % Save data for all detection times in continuous data
            for yy = 1:length(locs)
                
                tbeg = locs(yy) + tshiftc - win; % CCpoint + tp - window
                tbeg(tbeg<=0)=1;
                
                for vv = 1:length(tshiftc)
                    tend(vv) = locs(yy) + tshiftc(vv) + win;
                    if tend(vv) > size(datmpx,2); tend(vv) = size(datmpx,2); end;
                    
                    % In cells because of code legacy, could be now included in a classical array
                    datx{vv,1}=datmpx(vv,tbeg(vv):tend(vv));
                    daty{vv,1}=datmpy(vv,tbeg(vv):tend(vv));
                    datz{vv,1}=datmpz(vv,tbeg(vv):tend(vv));
                end
                
                tbegm = min(tbeg); tendm = max(tend);
                
                datxsv = dataxc(:,tbegm:tendm); % Get data from all stations just in case
                datysv = datayc(:,tbegm:tendm);
                datzsv = datazc(:,tbegm:tendm);
                
                hdr.st = 1/fsmaster; % Sample interval sec.
                yr = hd(71);
                datetmp = datenum([yr,1,1,0,0,0]);
                datetmp = datetmp + (hd(72)-1) + (hd(73)/24) +...
                    (hd(74)/1440) + ((hd(75)+(hd(76)/1000))/86400);
                
                hdr.tbeg = datetmp + (((tbeg-1)/fsmaster)/86400); % Beg. time of child event
                hdr.tbegsv = datetmp + (((tbegm-1)/fsmaster)/86400); % Beg. time of saved data
                hdr.tpk = peaks(yy); % Value of mad peak
                hdr.madthrs = madthres*xmad; % Threshold used for detection
                hdr.mad = madthres; % Value of mad level for detection
                hdr.compm = compsv; % Comps used for mad thres detection
                hdr.parentname = infopar(pp).name; % Name of parent events
                hdr.tbegcont = datetmp; % Beginning time of cont. data file
                hdr.txc = datetmp + ((locs(yy)/fsmaster)/86400); % Time of mad peak in cont data
                hdr.filt = [fcl fch]; % Filtering parameters in Hz
                
                filename = datestr(hdr.tbegsv,30);
                save([savepathtmp '/' filename],'datx','daty','datz',...
                    'hdr','hdrpar','datxsv','datysv','datzsv');
                
                clear tbeg* tend* yr hdr datetmp datx daty datz filename 
                clear datxsv datysv datzsv vv comp ind ind2
            end
        end
        
        % Continuous file # / parent file # / number of child events
        tracability((ii-1)*size(infopar,1)+pp,:) = [ii pp length(locs)];
%         disp(['Parent event ' dirname(1:end-4) ', # child events: ' ...
%             num2str(length(locs)) '. Par event ' num2str(pp) '/' num2str(size(infopar,1))])
        clear yy peaks locs datmpx datmpy datmpz ia ib lettre iv 
        clear dataxp datayp datazp xsum xmad hdrpar xcnew tshift*
        clear dataxs datays datazs parent savepathtmp dirname t0 var*
        clear stapc stasc kp ks tp ts stasv stasvp stasvs comp* stab
    end
    
    fprintf(fid,'%s\n',['Continuous file # ' num2str(ii) '/' num2str(size(infoc,1))]);
    disp(['Continuous file # ' num2str(ii) '/' num2str(size(infoc,1))])
    clear datax datay dataz hd kk sta jj pp dataxc datayc datazc stasv
    toc
end

% Save parameters used for this round of match filter in the same folder as
% where are located the folders of child events
save([savepath '/MatchFilterParams'])
fclose(fid);
