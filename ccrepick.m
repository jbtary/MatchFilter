% Repick events based on cross-correlation with data from parent events
% Child events are put at the location of the master event as a starting
% location. DT times are obtained using cross-correlation, then removing
% travel-times of the master event (tpicks - T0) and then finding the new
% T0 of the child event and updating picks with this new T0.
% 
% Cross-correlation is made using P- and S- wave picks of the master event
% to isolate these wavetrains and finding their max CC within child events.
% CC and DT are saved and used only if they are > to a minimum CC.
% 
% Child events with too little good CC differential times are moved to another
% folder. Events with too little number of connections will be removed as
% part of the hypoDD procedure, not in this program.
% 
% We are not using calculated travel-times because of the large misfits we
% get with the 1D velocity model (especially for far away stations).
% 
% To do:
% Use station corrections instead of shifting elevation of stations
% 
% INPUT:
% dirchild: path to folder containing the folders containing the child
% events (ex: '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ChildEvents/')
% dirpar: path to folder containing the matlab data files of the parent
% events (ex: '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ParDir/Picked/')
% newdir: new folder where to move events with # of CC data < minsta (ex:
% 'Infminsta')
% fs: data sampling frequency (same for everybody) (ex: 125 Hz)
% fcl: lower limit of pass-band filter (ex: 10 Hz)
% fch: higher limit of pass-band filter (ex: 30 Hz)
% winP: window of the P-wave to include in samples (ex: 125)
% winS: window of the S-wave to include in samples (ex: 190)
% mincc: min CC coef. to save it as a new pick for a given child event
% (between 0.5-1)
% minsta: min number of stations with mincc to keep the event in the
% reloction (>=3 would be recommended, but this selection can also be done
% by the hypoDD code)
% 
% OUTPUT:
% None

function ccrepick(dirchild,dirpar,newdir,mincc,minsta,fs,fcl,fch,winP,winS)

% Base selection on:
% - Cross-correlation coefficient > mincc
% - Number of stations with cc > mincc (minsta)

pardir = dir([dirchild 'Ev*']); % Get all folders of parent events

% First select events that are well correlated to the parent event
% Get their absolute times and DD time at the same time

for ii = 1:length(pardir) % Loop on parent events
    
    disp(['Parent # ' num2str(ii) ' ' pardir(ii).name])
    evlistmp = dir([dirchild pardir(ii).name '/*.mat']); % Get child events
    parent = load([dirpar pardir(ii).name '.mat']); % Get parent event
    
    % Data picked for parent event
    datparx = parent.dataxs;
    datpary = parent.datays;
    datparz = parent.datazs;
    compar = parent.comp;
    
    % T0 of parent event
    if isfield(parent,'newhdr') == 1 % In case the T0 has been corrected
        part0 = datenum(parent.newhdr(9),parent.newhdr(10),parent.newhdr(11),...
            parent.newhdr(12),parent.newhdr(13),parent.newhdr(14));
    else
        part0 = datenum(parent.hdr(9),parent.hdr(10),parent.hdr(11),parent.hdr(12),...
            parent.hdr(13),parent.hdr(14));
    end
    
    % Check if both fields exists
    if isfield(parent,'ppicks') == 0 && isfield(parent,'spicks') == 0
        disp(['No picks for ' pardir(ii).name '. Going to next event.'])
        continue;
    end
    
    % Check if both fields are empty
    if isempty(parent.ppicks) == 1 && isempty(parent.spicks) == 1
        disp(['No picks for ' pardir(ii).name '. Going to next event.'])
        continue;
    end
    
    for jj = 1:length(evlistmp)
        
        child = load([dirchild pardir(ii).name '/' evlistmp(jj).name]); % Get child event
        % Select matching data parent - child events (extracted data)
        compchild = child.hdr.compm;
        [~,ia,ib] = intersect(compchild,compar);
        
        tmpdatpx = datparx(ib,:); tmpdatpy = datpary(ib,:); tmpdatpz = datparz(ib,:);
        tmpdatcx = child.datx(ia,:); tmpdatcy = child.daty(ia,:); tmpdatcz = child.datz(ia,:); 
        
        % Get data with picks, separately for P and S waves
        compcomb = compar(ib); % Components with data in parent and child
        [~,iap,ibp] = intersect(compcomb,parent.stap);
        [~,ias,ibs] = intersect(compcomb,parent.stas);
        
        % Declare variables to avoid errors later
        ccp=[]; ccs=[]; tpc=[]; tsc=[]; staddp=[]; stadds=[]; 
        ccpt0=[]; ccst0=[];
        
        % P-waves part
        if isempty(iap) ~= 1
        Pdatpx = tmpdatpx(iap,:); % Data parent with P picks
        Pdatpy = tmpdatpy(iap,:); Pdatpz = tmpdatpz(iap,:);
        
        Pdatcx = tmpdatcx(iap,:); % Data child, corresponding components
        Pdatcy = tmpdatcy(iap,:); Pdatcz = tmpdatcz(iap,:);
        % Get corresponding P picks and corresponding stations
        ppicks = parent.ppicks(ibp,:); % Components re-ordered to match data automatically
        staddp = parent.stap(ibp,:);
        
        databeg = parent.datebeg(ib); databeg = databeg(iap);        
        ppicks = ((part0 + (ppicks/86400)) - databeg)*86400*fs; % Timing of picks in parent data in samples
        childbeg = child.hdr.tbeg(ia); childbeg = childbeg(iap);
        
        [ccp,tpc] = getddinfo(Pdatpx,Pdatpy,Pdatpz,Pdatcx,Pdatcy,Pdatcz,...
            ppicks,winP,childbeg,fcl,fch,fs); % Max CC coeff and abs times of picks
        
        % Get T0 for each station from new CC picks
        ccpt0 = tpc - (parent.ppicks(ibp,:)/86400);
        
        clear Pdatp* Pdatc* databeg childbeg
        end
        
        % S-waves part
        if isempty(ias) ~= 1
        Sdatpx = tmpdatpx(ias,:); % Data parent with S picks
        Sdatpy = tmpdatpy(ias,:); Sdatpz = tmpdatpz(ias,:);
        
        Sdatcx = tmpdatcx(ias,:); % Data child, corresponding components
        Sdatcy = tmpdatcy(ias,:); Sdatcz = tmpdatcz(ias,:);
        % Get corresponding S picks and corresponding stations
        spicks = parent.spicks(ibs,:);
        stadds = parent.stas(ibs,:);
        
        databeg = parent.datebeg(ib); databeg = databeg(ias);
        spicks = ((part0 + (spicks/86400)) - databeg)*86400*fs; % Timing of picks in parent data in samples
        childbeg = child.hdr.tbeg(ia); childbeg = childbeg(ias);
        
        [ccs,tsc] = getddinfo(Sdatpx,Sdatpy,Sdatpz,Sdatcx,Sdatcy,Sdatcz,...
            spicks,winS,childbeg,fcl,fch,fs);
        
        % Get T0 for each station from new CC picks
        ccst0 = tsc - (parent.spicks(ibs,:)/86400);
        
        clear Sdatp* Sdatc* databeg childbeg
        end
        
        ccmax = [ccp;ccs];
        cct0 = [ccpt0;ccst0];
        % Find new child T0 and update CC picks
        iv = find(ccmax > mincc);
        t0new = median(cct0(iv,1));
        
        hdr = child.hdr;
        hdr.t0child = t0new;
        
        hdr.ddppicks=[]; hdr.ddstap=[]; hdr.weightsp=[]; phaP=[];
        if isempty(iap) ~= 1
            ivp = find(ccp > mincc);
            hdr.ddppicks = (tpc(ivp,1) - t0new)*86400; % P-wave picks from CC to save
            hdr.ddstap = staddp(ivp,:); % Stations
            phaP = repmat('P',[length(hdr.ddppicks) 1]); % Phases
            hdr.weightsp = ccp(ivp,1); % Weights using CC coeffs
        end
        
        hdr.ddspicks=[]; hdr.ddstas=[]; hdr.weightss=[]; phaS=[];
        if isempty(ias) ~= 1
            ivs = find(ccs > mincc);
            hdr.ddspicks = (tsc(ivs,1) - t0new)*86400; % P-wave picks from CC to save
            hdr.ddstas = stadds(ivs,:); % Stations
            phaS = repmat('S',[length(hdr.ddspicks) 1]); % Phases
            hdr.weightss = ccs(ivs,1); % Weights using CC coeffs
        end
        clear ivp ivs ccmax cct0
        
        % Check if Ppicks are before Spicks (negative picks are removed
        % later in cctohypoDD
        kt = 1;
        for ll = 1:length(hdr.ddspicks)
            % Indice of possible match in hdr.ddstap
            test = strmatch(hdr.ddstas(ll,:),hdr.ddstap,'exact');
            spic = hdr.ddspicks(ll);
            ppic = hdr.ddppicks(test);
            if spic <= ppic % gather bad indexes
                ind(kt) = ll; kt=kt+1;
            end
            clear spic ppic test
        end
        
        if exist('ind','var')==1
            hdr.ddspicks(ind) = [];
            hdr.ddstas(ind,:) = [];
            phaS(ind) = [];
            hdr.weightss(ind) = [];
        end
        clear ind kt ll
        
        if length(iv) >= minsta % Save picks, T0 etc for current child event
            
            % New child header with new T0
            hdrtmp = child.hdrpar.hdr;
            % # YR MO DY HR MI SC LAT LON DEP MAG EH EZ RMS ID
            line = ['# ' datestr(t0new,'yyyy') ' ' datestr(t0new,'mm') ' ' ...
                datestr(t0new,'dd') ' ' datestr(t0new,'HH') ' ' datestr(t0new,'MM') ...
                ' ' datestr(t0new,'SS.FFF') ' ' num2str(hdrtmp(1,2),'%3.4f') ' ' ...
                num2str(hdrtmp(1,3),'%3.4f') ' ' num2str(hdrtmp(1,4),'%3.1f') ' ' ...
                num2str(hdrtmp(1,15),'%3.1f') ' ' num2str(10,'%3.1f') ' ' ...
                num2str(10,'%3.1f') ' ' num2str(hdrtmp(1,5),'%3.4f') ' ' ...
                num2str([evlistmp(jj).name(5:8) evlistmp(jj).name(10:15)],'%d')];
            
            hdrDD{1,1} = line;
            stadd = [hdr.ddstap;hdr.ddstas];
            picdd = [hdr.ddppicks;hdr.ddspicks];
            weidd = [hdr.weightsp;hdr.weightss];
            phadd = [phaP;phaS];
            
            % Sorting data by station number and Ifremer/Koeri for hypoDD
            [~,ind2] = sort(str2num(stadd(:,2:3)),1); % Station number
            [~,ind4] = sort(stadd(ind2,1),1,'descend'); % K and G sorting
            ind5 = ind2(ind4);
            
            stadd = stadd(ind5,:); picdd = picdd(ind5,:);
            weidd = weidd(ind5,:); phadd = phadd(ind5,:);
            
            for nn = 1:length(picdd)
                statmp = stadd(nn,:);
                if strcmp(statmp(1),'G')==1; statmp=['OB' statmp(2:3)]; end;
                if strcmp(statmp(1),'K')==1; statmp=['KG' statmp(2:3)]; end;
                
                hdrDD{nn+1,1} = [statmp ' ' num2str(picdd(nn,1),'%3.6f') ...
                    ' ' num2str(weidd(nn,1),'%1.1f') ' ' phadd(nn,1)];
                
                clear statmp
            end
            hdr.DDdata = hdrDD;
            
            save([dirchild pardir(ii).name '/' evlistmp(jj).name],'hdr','-append');
            
            clear ind* stadd picdd weidd phadd nn hdr hdrDD hdrtmp line
            
        else
            % Move event to newdir folder
            namedir = [dirchild pardir(ii).name '/' newdir];
            if exist(namedir,'dir') ~= 7
                mkdir(namedir);
            end
            movefile([dirchild pardir(ii).name '/' evlistmp(jj).name],...
                namedir);
            
            clear namedir
        end
        
        clear tmpdat* compchild ia ib child ttimes t0diff t0new kk
        clear ccp* ccs* compcomb iap ias ibp ibs iv phaP phaS ppicks spicks
        clear staddp stadds tpc tsc
    end
    
    clear datpar* parent evlistmp compar part0
end

save([dirchild '/ccrepick_params'])

% Get diff times and cc coeffs for any input data
function [ccm,tps] = getddinfo(datpx,datpy,datpz,datcx,datcy,...
    datcz,picks,win2,cbeg,fcl,fch,fs)
% picks in samples
% win2 in samples

fNy = fs/2;
[b,a]=butter(6,[fcl fch]/fNy,'bandpass');

for kk = 1:length(datcx) % Loop on the stations
    
    % Filter data and normalize
    datcx2 = detrend(datcx{kk,:}'); % Child event data
    datcx2 = filtfilt(b,a,datcx2); datcx2 = datcx2';
    datcy2 = detrend(datcy{kk,:}');
    datcy2 = filtfilt(b,a,datcy2); datcy2 = datcy2';
    datcz2 = detrend(datcz{kk,:}');
    datcz2 = filtfilt(b,a,datcz2); datcz2 = datcz2';
    
    datcx2=datcx2/max(abs(datcx2));
    datcy2=datcy2/max(abs(datcy2));
    datcz2=datcz2/max(abs(datcz2));
    
    datpx2 = detrend(datpx{kk,:}'); % Parent event data
    datpx2 = filtfilt(b,a,datpx2); datpx2 = datpx2';
    datpy2 = detrend(datpy{kk,:}');
    datpy2 = filtfilt(b,a,datpy2); datpy2 = datpy2';
    datpz2 = detrend(datpz{kk,:}');
    datpz2 = filtfilt(b,a,datpz2); datpz2 = datpz2';
    
    datpx2=datpx2/max(abs(datpx2));
    datpy2=datpy2/max(abs(datpy2));
    datpz2=datpz2/max(abs(datpz2));
    
    % Timing of arrival
    t1 = round(picks(kk)-(win2/2));
    t2 = round(picks(kk)+win2);
    t1(t1<=0) = 1; t2(t2>length(datpx2)) = length(datpx2);
    % Avoid errors with bad picks outside the data (negative t2)
    if t2 < 5; tps(kk,1) = NaN; ccm(kk,1) = NaN; continue; end;
    pdatx = datpx2(1,t1:t2); % Selection of arrival
    pdaty = datpy2(1,t1:t2);
    pdatz = datpz2(1,t1:t2);
    
    % Cross-correlation to estimate coeff and timing
    % Next is searching for parent in child data
    [xcx,~] = xcorr(datcx2,pdatx);
    [xcy,~] = xcorr(datcy2,pdaty);
    [xcz,lag] = xcorr(datcz2,pdatz);
    
    xc3 = [xcx;xcy;xcz];
    [iv,~] = find(abs(xc3) == max(max(abs(xc3))));
    if length(iv) > 1; iv = iv(1); end;
    % Refine DD time by quadratic interpolation
    xctmp = xc3(iv,:);
    
    [~,lc1]=max(abs(xctmp));
    if lc1==1; lc1=2; end; % Avoid problem for next line
    if lc1==length(xctmp); lc1=length(xctmp)-1; end;
    P=polyfit([lag(lc1-1) lag(lc1) lag(lc1+1)],[xctmp(lc1-1) xctmp(lc1) xctmp(lc1+1)],2);
    dt = (-P(2)/(2*P(1)))/fs; clear P
    clear xcx xcy xcz xc3
    
    tps(kk,1) = cbeg(kk) + (dt/86400); % Absolute time of CC max
    
    % Use common part with parent event to calculate cross-corr
    % coefficient (no timing here because there might be little
    % differences with the real timing of the events)
    % CHECK if timing, selection and size are corrects (yes)
    s1 = lc1-5-length(datcx2); s1(s1<=0) = 1; 
    s2 = lc1-5-length(datcx2)+length(pdatx)-1; 
    s2(s2>length(datcx2))=length(datcx2);
    [xcx,~] = xcorr(datcx2(1,s1:s2),pdatx(1,1:s2-s1+1),'coeff');
    [xcy,~] = xcorr(datcy2(1,s1:s2),pdaty(1,1:s2-s1+1),'coeff');
    [xcz,~] = xcorr(datcz2(1,s1:s2),pdatz(1,1:s2-s1+1),'coeff');
    clear s1 s2
    
    % Get maximum CC value on same component as timing
    xc3 = [xcx;xcy;xcz];
    ccm(kk,1) = max(abs(xc3(iv,:)));
    
    clear datpx2 datpy2 datpz2 xc* lc1 lag
    clear datcx2 datcy2 datcz2 iv pdat* t1 t2 dt
end

clear fNy b a ii


% %%%%% To move files up a directory from newdir
% for ii = 1:length(pardir)
% movefile([dirchild pardir(ii).name '/' newdir '/*.mat'],[dirchild pardir(ii).name '/']);
% end
% 
% %%%%% To count the number of events remaining
% pardir = dir([dirchild 'Ev*']);
% evlist = []; dirlist = [];
% for ii = 1:length(pardir)
%     evlistmp = dir([dirchild pardir(ii).name '/*.mat']);
%     evlist = [evlist;evlistmp];
%     dirlistmp = repmat(pardir(ii).name,[length(evlistmp) 1]);
%     dirlist = [dirlist;dirlistmp];
%     clear evlistmp dirlistmp
% end
