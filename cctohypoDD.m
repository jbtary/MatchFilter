% Format parent events and child events picks and header into an input
% hypoDD file
% 
% Input:
% dirchild: path to folder containing the folders containing the child
% events (ex: '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ChildEvents/')
% dirpar: path to folder containing the matlab data files of the parent
% events (ex:
% '/media/jbtary/DeLaFeuille/BackUP/Marmara/Data/ParDir/Picked/'), not in
% use for now
% filename: name of the text file where the input data will be saved. This
% file will be saved in the dirchild directory.
% idcount: number to put in front of IDs of some events when there are some
% repeating IDs (impossible in HypoDD). The idcount variable will decrease
% as some repeating IDs are found. (ex: 999) Keep in mind that IDs in
% HypoDD can't be too long either (~10 char max)
% 
% Output:
% ddfile: cell array containing data formatted for hypoDD
% keeptrack: number of child events per parent events and number of IDs
% duplicates
% Text file for hypoDD input

function [ddfile,keeptrack] = cctohypoDD(dirchild,filename,idcount)

pardir = dir([dirchild 'Ev*']); % Get all folders of parent events

idcsv = idcount; % Save number of duplicate IDs possible
count = 1; % ids = [];
for ii = 1:length(pardir)
    
    evlistmp = dir([dirchild pardir(ii).name '/*.mat']); % Get child events
    disp(['Parent # ' num2str(ii) ' ' pardir(ii).name ', # of child events: ' ...
        num2str(length(evlistmp))])
    keeptrack{ii,1} = ['Parent # ' num2str(ii) ' ' pardir(ii).name ', # of child events: ' ...
        num2str(length(evlistmp))];
    
    if isempty(evlistmp) == 1 % If no child events, go to next parent event
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parent data part of code is going here %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj = 1:length(evlistmp)
        
        child = load([dirchild pardir(ii).name '/' evlistmp(jj).name],'hdr'); % Get child event
        % Header line of child event
        ddfile{count,1} = child.hdr.DDdata{1,1};
        count = count + 1;
        
        % Phase data of child event (remove negative picks)
        for kk = 2:length(child.hdr.DDdata)
            linetmp = child.hdr.DDdata{kk,1};
            nbs = str2num(linetmp(6:end-1));
            
            % Only save positive picks (ideally, P picks should also be < S
            % picks and take into account CC of P to know if it's a good pick)
            if sign(nbs(1)) > 0
                ddfile{count,1} = linetmp;
                count = count+1;
            end
            
            clear linetmp nbs
        end
        
        clear child
    end
    
    % In case there are lots of duplicate IDs we reinitialize the idcount when
    % it reaches 0. It's likely that reinitializing will not create
    % duplicate IDs, but it's a possibility.
    if idcount == 0
       idcount = idcsv;
       disp('idcount reinitialized, lots of duplicate IDs')
    end
    
    clear evlistmp parent jj kk
end
disp(['Number of duplicates: ' num2str(idcsv - idcount)])
clear count
keeptrack{ii+1,1} = ['Number of duplicates: ' num2str(idcsv - idcount)];

% Create input file
fid = fopen(filename,'w'); % new file
for ii = 1:length(ddfile)
    fprintf(fid,'%s\n',ddfile{ii,1});
end
fclose(fid);

% Move file to final destination
movefile(filename,dirchild);

clear idsv

%%%%%%%%%%%%%%%%%%%
% Parent files part
% parent = load([dirpar pardir(ii).name '.mat'],'newhdr','hDDph'); % Parent data
% 
% yy = num2str(parent.newhdr(9)); mm = num2str(parent.newhdr(10),'%02d');
% dd = num2str(parent.newhdr(11),'%02d');
% HH = num2str(parent.newhdr(12),'%02d'); MM = num2str(parent.newhdr(13),'%02d');
% SS = num2str(parent.newhdr(14),'%02.3f');
% % Avoid EH and EZ in meters instead of kms
% if parent.newhdr(6) > 500 || parent.newhdr(7) > 500 || parent.newhdr(8) > 500
%     parent.newhdr(6) = parent.newhdr(6)/1000;
%     parent.newhdr(7) = parent.newhdr(7)/1000;
%     parent.newhdr(8) = parent.newhdr(8)/1000;
% end
% EH = num2str(sqrt(parent.newhdr(6)^2 + parent.newhdr(7)^2),'%3.1f');
% EZ = num2str(parent.newhdr(8),'%3.1f');
% 
% % Avoid problem with repeating IDs (only for parent events)
% iv = find(ids == parent.newhdr(1));
% if isempty(iv) ~= 1
%     id = [num2str(idcount) num2str(parent.newhdr(1),'%g')];
%     idcount = idcount -1;
% else
%     id = num2str(parent.newhdr(1),'%g');
% end
% ids(ii) = parent.newhdr(1); % Inventory of parent IDs
% 
% line = ['# ' yy ' ' mm ' ' dd ' ' HH ' ' MM ' ' SS ' ' ...
%     num2str(parent.newhdr(2),'%3.4f') ' ' num2str(parent.newhdr(3),'%3.4f') ' ' ...
%     num2str(parent.newhdr(4),'%3.1f') ' ' num2str(parent.newhdr(15),'%1.1f') ' ' ...
%     EH ' ' EZ ' ' num2str(parent.newhdr(5),'%2.3f') ' ' id];
% % Header line for parent event
% ddfile{count,1} = line;
% clear yy mm dd HH MM SS EH EZ line id idtmp
% 
% % Add phase data from parent event
% count = count+1;
% for jj=1:length(parent.hDDph)
%     linetmp = parent.hDDph{jj,1};
%     if strcmp(linetmp(1),'G')==1; linetmp=['OB' linetmp(2:end)]; end;
%     if strcmp(linetmp(1),'K')==1; linetmp=['KG' linetmp(2:end)]; end;
%     
%     ddfile{count,1} = linetmp;
%     count = count+1;
%     clear linetmp
% end
