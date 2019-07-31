

function img=read_UFXC_ubit32_ParSpars(varargin)

file_name = varargin{1};
kinetics_mode_yes_no = varargin{2};

fr_leng_CW=11;
fr_leng_Burst=11;

tic

m = memmapfile(file_name, 'Offset',0,'Repeat',Inf,'Writable',false,'Format','uint32');
data = m.Data;
whos data;

disp('Begin Picking Bits');

if (kinetics_mode_yes_no == 1)

    pix_burst_batch=double(bitshift(data,fr_leng_Burst-32));  % Piotr said it shouldn't start from 0    
    pos=find((pix_burst_batch(2:end)-pix_burst_batch(1:end-1))<-2);
    
    for i=1:numel(pos)
        pix_burst_batch(pos(i)+1:end) = pix_burst_batch(pos(i)+1:end)+2^fr_leng_Burst;    
    end
    
    
%%   
    
    pix_burst_batch=pix_burst_batch-pix_burst_batch(1);
      
    pix_burst_batch=bitshift(pix_burst_batch,-1);
%%
%     
    pix_burst_frame=double(bitand(hex2dec('F'),bitshift(data,-17)))-1;
    
    pix_frame=pix_burst_frame+12*pix_burst_batch;
    pix_count=double(bitand(3,bitshift(data,-15)));
    pix_ind=double(bitand(hex2dec('7FFF'),data))+1;
    
    img=sparse(pix_ind,pix_frame,pix_count,128*256,12*(pix_burst_batch(end)+1));
    
else

    pix_frame_orig=double(bitshift(data,fr_leng_CW-32))+1;    
    pix_count=double(bitand(3,bitshift(data,-15)));
    pix_ind=double(bitand(hex2dec('7FFF'),data))+1;
 

%% Qingteng, Suresh and Piotr added on 08/14, Qingteng added comments on 10/09

    % This approach fixes scrambles of frame index close to 2048 
    % (e.g. ... 2047 2047 2048 2048 1 1 2048 1 2048 1 1 2 2 2 3 3 ...)
    % caused by time differences in CPUs in arriving at the same checkpoint

    % As long as block 100, 200, 300 and 1500 don't mix (like frame 232
    % from one CPU doens't get in between frame 300~1500 of the other CPUs), 
    % this approach should always work.
    
    dim=numel(pix_frame_orig);

    BF=0;   % Multiples of 2048 for incrementing frame numbers
    SM=0;   % Marker for determining when to increment frames

    time_stamp=zeros(dim,1);

    for ii=1:dim

        if SM == 0

            time_stamp(ii)=pix_frame_orig(ii)+BF;  % Increment using BF when SM = 0.
            if pix_frame_orig(ii)>1500  % SM flipped to 1 beyond pix_frame_orig reaching 1500 (Max=2048)
                SM=1;
            end
            
        elseif SM == 1  % After SM flipped to 1, enter this block

            if pix_frame_orig(ii)<300
                time_stamp(ii)=pix_frame_orig(ii)+BF+2048;  % For pix_frame_orig beyond 1500 but hasn't reached 300, increment 2048 with this block
            else
                time_stamp(ii)=pix_frame_orig(ii)+BF; % For pix_frame_orig beyond 300, increment using BF. BF is updated in the next parallel if block
            end

            if pix_frame_orig(ii)>100 && pix_frame_orig(ii)< 200 % For pix_frame_orig between 100 and 200, increment BF by 2048, set SM back to 0.
                SM=0;
                BF=BF+2048;
            end

        end

    end
  
%% Original

%     [~,ind]=sort(time_stamp);    
%     pix_frame=time_stamp-time_stamp(1)+1;    
%     pix_ind=pix_ind(ind);
%     pix_count=pix_count(ind);

%% Qingteng changed on 08/19 2:23 am    

    [pix_frame,ind]=sort(time_stamp);
    pix_frame=pix_frame-pix_frame(1)+1;
    pix_ind=pix_ind(ind);
    pix_count=pix_count(ind);
        

%%
    
        
    img=sparse(pix_ind,pix_frame,pix_count,128*256,max(pix_frame));

end



disp('File reading is done');

toc

% 
%  pix=bitand(hex2dec('7FFF'),udata)+1;
%  cts=bitand(3,bitshift(data,-15));
%  
%   fr=bitshift(data,-21);

