% Remove duplicates
% Detection based on cross-correlation for all master events is ending up
% with many duplicates due to some master events being highly similar. This
% has the advantage to spanning more waveforms but the disadvantage of
% generating duplicates. There are also duplicates of the master events
% themselves, when some are also highly correlated or the master event
% itself.
% Removal of duplicates based on timing of events and then number of
% stations and MAD threshold/peak.
% 
% Input:
% dirchild - directory with folders of master events containing their child
% events (ex: '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ChildEvents/')
% dupldir - directory name to move duplicates into (inside each parent
% event folder, ex: 'Duplicates')
% minsamps - minimum number of seconds between events (ex: =ntime/fsmaster
% in MatchFilt_Marmara.m = 10)
% 
% Output:
% finallist - list of events kept (non-duplicates), associated parent
% event, associated duplicates, associated parent events of duplicates

function finallist = rm_duplicates(dirchild,dupldir,minsamps)

pardir = dir([dirchild 'Ev*']); % Get all folders of parent events

% Generate complete child event list and associated directories
evlist = []; dirlist = [];
for ii = 1:length(pardir)
    evlistmp = dir([dirchild pardir(ii).name '/*.mat']);
    evlist = [evlist;evlistmp];
    dirlistmp = repmat(pardir(ii).name,[length(evlistmp) 1]);
    dirlist = [dirlist;dirlistmp];
    clear evlistmp dirlistmp
end

dates = arrayfun(@(x) x.name(1:15),evlist,'UniformOutput',false);
dates = char(dates);
evlistnames = char(extractfield(evlist,'name')');
datenb = datenum(dates,'yyyymmddTHHMMSS'); % Beginning time of events from file names

% Remove events that are at the time of any master event
datemaster = arrayfun(@(x) x.name(10:end),pardir,'UniformOutput',false);
datemaster = datenum(char(datemaster),'yymmdd_HHMMSS');
for ii = 1:length(datemaster) % Loop on master events with child events
    iv = find(datenb>(datemaster(ii)-(minsamps/86400)) & datenb<(datemaster(ii)+(minsamps/86400)));
    
    if isempty(iv) ~= 1 % If there are some duplicates
        for jj = 1:length(iv) % Loop on duplicates
            namedir = [dirchild dirlist(iv(jj),:) '/' dupldir];
            if exist(namedir,'dir') ~= 7
                mkdir(namedir);
            end
            movefile([dirchild dirlist(iv(jj),:) '/' evlistnames(iv(jj),:)],...
                namedir);
                        
            clear namedir
        end
        % Remove duplicates from lists
        datenb(iv) = [];
        dates(iv,:) = [];
        evlist(iv) = [];
        dirlist(iv,:) = [];
        evlistnames(iv,:) = [];
    end
    clear iv jj
end

[h,~,bins] = histcounts(datenb,round((max(datenb - min(datenb)))/(minsamps/86400)));
[~,indh] = find(h > 1);

for ii = 1:length(indh) % Loop on number of events with duplicates
    
    nametmp = find(bins == indh(ii));
    for jj = 1:length(nametmp) % Get hdr info from each event in this duplicates group
        load([dirchild dirlist(nametmp(jj),:) '/' evlist(nametmp(jj)).name],'hdr')
        % Number of stations used for MAD, MAD value, MAD thres, difference
        % MAD value-MAD thres
        evlistmp(jj,:) = [length(hdr.compm) hdr.madthrs hdr.tpk hdr.tpk-hdr.madthrs];
        
        clear hdr
    end
    
    % Sort events depending first on the number of stations and then on the
    % difference between the MAD threshold and the observed MAD value
    [~,indsrt] = sortrows(evlistmp,[4 1]);
    
    % Save names of other potential parent events in hdr, as well as names
    % of duplicates
    filesv = evlist(nametmp(indsrt(1))).name;
    load([dirchild dirlist(nametmp(indsrt(1)),:) '/' filesv],'hdr')
    hdr.evdupl = evlistnames(nametmp(indsrt(2:end)),:);
    hdr.pardupl = dirlist(nametmp(indsrt(2:end)),:);
    save([dirchild dirlist(nametmp(indsrt(1)),:) '/' filesv],'hdr','-append');
    
    % Move duplicates to another directory
    for jj = 1:length(indsrt)-1
        namedir = [dirchild dirlist(nametmp(indsrt(jj+1)),:) '/' dupldir];
        if exist(namedir,'dir') ~= 7
            mkdir(namedir);
        end
        movefile([dirchild dirlist(nametmp(indsrt(jj+1)),:) '/' evlistnames(nametmp(indsrt(jj+1)),:)],...
            namedir);
        
        clear namedir
    end
    
    finallist{ii,1} = dirlist(nametmp(indsrt(1)),:); % Parent event of saved event
    finallist{ii,2} = filesv; % Name of saved event
    finallist{ii,3} = hdr.pardupl; % Parent events of duplicate events
    finallist{ii,4} = hdr.evdupl; % Name of duplicate events
    
    clear hdr filesv indsrt jj nametmp evlistmp indsrt
end

save([dirchild '/ListDuplicates'])
