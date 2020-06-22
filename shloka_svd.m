clc
clear all

%%
pure_swaras = ['0905'; '0906'; '0907'; '0908'; '0909' ; '090A';'090B';'090C';'090F';'0910'...
    ;'0913';'0914'];
matras = ['093D';'093E';'093F';'0940';'0941';'0942';'0943';'0944';'0947';'0948';'094B';'094C'];
am_aha =['0902';'0903'];
v1 = ['0915'; '0916';'0917';'0918';'0919'];
v2 = ['091A';'091B';'091C';'091D';'091E'];
v3 = ['091F';'0920';'0921';'0922';'0923'];
v4 = ['0924';'0925';'0926';'0927';'0928'];
v5 = ['092A';'092B';'092C';'092D';'092E'];
av = ['092F';'0930';'0932';'0933';'0935';'0936';'0937';'0938';'0939'];
v =[v1;v2;v3;v4;v5;av];
halanta = ['094D'];
local_default = ['wxyz';'wxyz';'wxyz';'wxyz';'wxyz'];
%%

%fid=fopen('devanagari.txt','rt');
fid=fopen('vsn_devn.txt','rt');

bigspace1 =[]; bigspace2 = []; bigspace3 = [];

while ~feof(fid)
    Maze=[];
    Maze2=[];
    split_matrix=[];
    am_aha_matrix=[];
    Line=fgetl(fid);
    Maze=unicode2native(Line);
    Maze=dec2hex(Maze);
    
    size(Maze)
    if mod(size(Maze,1),2)~=0
        Maze=Maze(1:end-1,:);
    end
    %Maze = Maze(3:end-1,:); %to remove new line, cr and end 00
    %size(Maze)
    
    if isempty(Maze) % to get rid of empty lines in between actual lines
        disp('lineis empty');
        continue
    end
    
    Maze2 = [Maze(1:2:end,:) Maze(2:2:end,:)];
    local_index = 1;
    local = local_default;
    clear_local = 0;
    for i=1:length(Maze2)
        given = Maze2(i,:);
        
        
        if clear_local == 1
            local = local_default;
            clear_local = 0;
        end
        
        if ismember(given,pure_swaras,'rows')
            disp(['Pure swara' given])
            local(local_index,:) = given;
            local_index=1;
            
            next = Maze2(i+1,:);
            if ismember(next,am_aha,'rows')
                am_aha_matrix = [am_aha_matrix; next];
            else
                am_aha_matrix = [am_aha_matrix; 'wxyz'];
            end
            
            
        elseif ismember(given, halanta,'rows')
            disp('Is halanta')
            %local = ['1234';'1234';'1234';'1234';'1234'];
            %local_index=1; 
            continue
        elseif ismember(given,matras,'rows')
            disp('Is matra')
            local = local_default;
            local_index=1;
            next = Maze2(i+1,:);
            if ismember(next,am_aha,'rows')
                am_aha_matrix = [am_aha_matrix; next];
            else
                am_aha_matrix = [am_aha_matrix; 'wxyz'];
            end
            
            continue
        elseif ismember(given,v,'rows')
            disp('Is vyanjana')
            next = Maze2(i+1,:);
            if ~ismember(next,halanta,'rows')
                disp('Next is not halanta')
                if ismember(next, matras, 'rows')
                    disp('Next is matra')
                    local(local_index,:) = given;
                    local_index=local_index+1;
                    local(local_index,:)=next; % CHANGE MATRAS TO ACTUAL TEXT
                    %local_index=1;
                    
%                     i=i+1;
                else
                    disp('Next is not matra')
                    local(local_index,:) = given;
                    local_index=local_index+1;
                    local(local_index,:)='0905'; % CHANGE MATRAS TO ACTUAL TEXT
                    local_index=1;
                    clear_local=1;
                    
                    next2 = Maze2(i+1,:);
                    disp(['next2 = ' next2])
                    if ismember(next2,am_aha,'rows')
                        am_aha_matrix = [am_aha_matrix; next];
                    else
                        am_aha_matrix = [am_aha_matrix; 'wxyz'];
                    end
                    
                end
            else    % is halanta
                disp('Next is halanta')
                local(local_index,:) = given;
                local_index = local_index+1;
                continue
%                 i=i+1;
            end
        else
            disp(['not swara not vyanjana' given])
            continue
                    
                    
        end
        disp(reshape(local',[1,4*5]))
        split_matrix=[split_matrix ; reshape(local',[1,4*5])];
        while(split_matrix(end,end)=='z')
            split_matrix(end,:)=circshift(split_matrix(end,:),4);
        end
        
        %replacing mataras by pure swaras
        for j = 1:length(matras)
            if split_matrix(end,end-3:end) == matras(j,:)
                split_matrix(end,end-3:end) = pure_swaras(j,:);
            end
        end
    end

    %Now replace everything with their indices
    subspace1 =  zeros(length(split_matrix),4);     %vyanjanas
    subspace2 = zeros(length(split_matrix),1);        %swaras
    subspace3 = zeros(length(split_matrix),1);        %am_aha
    for j=1:4 %four vyanjanas, horizontal index
        for j2 = 1:length(split_matrix) %16 chars, vertical index
            for j3 = 1:length(v) %to iterate through all vyanjanas
                disp([split_matrix(j2, (j-1)*4 +1:(j-1)*4 +4) v(j3,:)])
                if prod(int8(split_matrix(j2, (j-1)*4 +1:(j-1)*4 +4) == v(j3,:)),'all')==1
                    subspace1(j2,j)= j3;
                    break
                end
            end
        end
    end
    
    for j = 1:length(subspace2) %swaras
        for j2 = 1:length(pure_swaras)
            if prod(split_matrix(j,end-3:end) == pure_swaras(j2,:),'all') == 1
                subspace2(j) = j2;
                break
            end
        end
    end
    
    for j = 1:length(subspace3) %am_aha
        for j2 = 1:length(am_aha(:,1))
            if prod(am_aha_matrix(j,:) == am_aha(j2, :), 'all')  == 1
                subspace3(j) = j2;
                break
            end
        end
    end
    
    bigspace1 = cat(3, bigspace1, subspace1);
    bigspace2 = [bigspace2 subspace2];
    bigspace3 = [bigspace3 subspace3];
    
    
end

%safety duplications
%bsp1 = bigspace1; bsp2 = bigspace2; bsp3 = bigspace3;

mean_BS1 = mean(bigspace1,3);
mean_BS2 = mean(bigspace2,2);
mean_BS3 = mean(bigspace3,2);

%should we normalize here? with the length of v? or normalize U with its
%own length?
bigspace1 = (bigspace1 - mean_BS1)./length(v(:,1));
bigspace2 = (bigspace2 - mean_BS2)./length(pure_swaras(:,1));
bigspace3 = (bigspace3 - mean_BS3)./length(am_aha(:,1));

%%
[U_BS11,S_BS11,V_BS11] = svd(squeeze(bigspace1(:,1,:)));
[U_BS12,S_BS12,V_BS12] = svd(squeeze(bigspace1(:,2,:)));
[U_BS13,S_BS13,V_BS13] = svd(squeeze(bigspace1(:,3,:)));
[U_BS14,S_BS14,V_BS14] = svd(squeeze(bigspace1(:,4,:)));

[U_BS2,S_BS2,V_BS2] = svd(bigspace2(:,1,:));
[U_BS3,S_BS3,V_BS3] = svd(bigspace3(:,1,:));

%%
%Principal Components
PC_BS11 = abs(round(U_BS11(:,1)*length(v(:,1)) + mean_BS1(32,1)));
PC_BS12 = abs(round(U_BS12(:,1)*length(v(:,1)) + mean_BS1(32,2)));
PC_BS13 = abs(round(U_BS13(:,1)*length(v(:,1)) + mean_BS1(32,3)));
PC_BS14 = abs(round(U_BS14(:,1)*length(v(:,1)) + mean_BS1(32,4)));

PC_BS2 = abs(round(U_BS2(:,1)*length(pure_swaras(:,1)) + mean_BS2));
PC_BS3 = abs(round(U_BS3(:,1)*length(am_aha(:,1)) + mean_BS3));

%limit their max values
PC_BS11(PC_BS11>length(v(:,1))) = length(v(:,1));
PC_BS12(PC_BS12>length(v(:,1))) = length(v(:,1));
PC_BS13(PC_BS13>length(v(:,1))) = length(v(:,1));
PC_BS14(PC_BS14>length(v(:,1))) = length(v(:,1));
PC_BS2(PC_BS2>length(pure_swaras(:,1))) = length(pure_swaras(:,1));
PC_BS3(PC_BS3>length(am_aha(:,1))) = length(am_aha(:,1));

%%
concise=[];
for i=1:length(PC_BS11)
    
    add_halanta = false;
    %vyanjanas
    if PC_BS11(i)~=0
        concise = [concise v(PC_BS11(i),:)];        
        add_halanta = true;
    end
    if PC_BS12(i)~=0
        if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(PC_BS12(i),:)];
        add_halanta = true;
    end
    if PC_BS13(i)~=0
         if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(PC_BS13(i),:)];
        add_halanta = true;
    end
    if PC_BS14(i)~=0
         if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(PC_BS14(i),:)];
        add_halanta = true;
    end
    
    %swaras
    if PC_BS2(i)~=0
        %if add_halanta is true, it means a vyanjana has preceeded the
        %swara. In that case, we replace the swara by its matra. Otherwise,
        %the swara remains.
        if add_halanta == true
            if PC_BS2(i) ~= 1 %should not add anything for a
                concise = [concise matras(PC_BS2(i),:)];
            %else - add nothing
            end
        else
            concise = [concise pure_swaras(PC_BS2(i),:)];
        end
    end
    
    %am_aha
    if PC_BS3(i)~=0 && PC_BS2(i)~=0 %do not add am_aha without swaras!
        concise = [concise am_aha(PC_BS3(i),:)];
    end
end
%%
concise_dec = [concise(1:2:end)' concise(2:2:end)'];
concise_dec = ['FE'; 'FF'; concise_dec;'00';'00'];
concise_dec = hex2dec(concise_dec);
concise_dec = [concise_dec(1:2:end) concise_dec(2:2:end)];

%%
fid =fopen('output.txt','w');
uni = []; %native2unicode(concise_dec);

for i = 1:length(concise_dec(:,1))
    uni = [uni; native2unicode(concise_dec(i,:))];
    concise_dec(i,:)
    %fwrite(fid,uni(i,:));
    fprintf(fid,'%0s%0s%0s%0s',uni(i,:));
end

fclose(fid);
%%
%now find the coeffs - U'(X-mu), mu = mean
%means are already subtracted from X in bigspaceN values
%of the form 32xN - 32 = no of syllables per line, N = number of lines
coeffs_BS11 = U_BS11' * squeeze(bigspace1(:,1,:)); 
coeffs_BS12 = U_BS12' * squeeze(bigspace1(:,2,:));
coeffs_BS13 = U_BS13' * squeeze(bigspace1(:,3,:));
coeffs_BS14 = U_BS14' * squeeze(bigspace1(:,4,:));
coeffs_BS2  = U_BS2'  * bigspace2;
coeffs_BS3  = U_BS3'  * bigspace3;

%reconstruction of first line for example
r = 32; %required rank
%every column of recon line corresponds to the split_marix equivalent for
%that line
recon_line = [sum(coeffs_BS11(1:r,1).* bigspace1(:,1,1)',1)' sum(coeffs_BS12(1:r,1).* bigspace1(:,2,1)',1)' sum(coeffs_BS13(1:r,1).* bigspace1(:,3,1)',1)' sum(coeffs_BS14(1:r,1).* bigspace1(:,4,1)',1)' sum(coeffs_BS2(1:r,1).* bigspace2(:,1)',1)' sum(coeffs_BS3(1:r,1).* bigspace3(:,1)',1)'];
%add the means now
recon_line = recon_line + [mean_BS1 mean_BS2 mean_BS3];
recon_line = round(recon_line)


%repeat file writing steps

concise=[];
for i=1:length(recon_line(:,1))
    
    add_halanta = false;
    %vyanjanas
    if recon_line(i,1)~=0
        concise = [concise v(recon_line(i,1),:)];        
        add_halanta = true;
    end
    if recon_line(i,2)~=0
        if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(recon_line(i,2),:)];
        add_halanta = true;
    end
    if recon_line(i,3)~=0
         if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(recon_line(i,3),:)];
        add_halanta = true;
    end
    if recon_line(i,4)~=0
         if add_halanta == true
            concise = [concise halanta];
        end
        concise = [concise v(recon_line(i,4),:)];
        add_halanta = true;
    end
    
    %swaras
    if recon_line(i,5)~=0
        %if add_halanta is true, it means a vyanjana has preceeded the
        %swara. In that case, we replace the swara by its matra. Otherwise,
        %the swara remains.
        if add_halanta == true
            if recon_line(i,5) ~= 1 %should not add anything for a
                concise = [concise matras(recon_line(i,5),:)];
            %else - add nothing
            end
        else
            concise = [concise pure_swaras(recon_line(i,5),:)];
        end
    end
    
    %am_aha
    if recon_line(i,6)~=0 && recon_line(i,5)~=0 %do not add am_aha without swaras!
        concise = [concise am_aha(recon_line(i,6),:)];
    end
end
%%
concise_dec = [concise(1:2:end)' concise(2:2:end)'];
concise_dec = ['FE'; 'FF'; concise_dec;'00';'00'];
concise_dec = hex2dec(concise_dec);
concise_dec = [concise_dec(1:2:end) concise_dec(2:2:end)];

%%
fid =fopen('output_recon.txt','w');
uni = []; %native2unicode(concise_dec);

for i = 1:length(concise_dec(:,1))
    uni = [uni; native2unicode(concise_dec(i,:))];
    concise_dec(i,:)
    %fwrite(fid,uni(i,:));
    fprintf(fid,'%0s%0s%0s%0s',uni(i,:));
end

fclose(fid);

