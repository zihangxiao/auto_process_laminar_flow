% read all csvs in the dir (The only thing to change)
exp_store_dir = '/Users/51325/Desktop/caltech_analysis/all_raw_data_position/10%w BE/';
out_dir = '/Users/51325/Desktop/caltech_analysis/all_raw_data_position/10%w BE/';

% find result files
exps = dir(fullfile(exp_store_dir,'*.csv'));
%load all files with specific document type
expnames = {exps.name};

%extract name
all_unique_exp = cell(1,size(expnames,2));
%size：return row/collum number，1 for row, 2 for collum

for k = 1 : size(expnames,2)
    % add exp name
    tmp = strsplit(expnames{1,k},'_');
    all_unique_exp{1,k} = strjoin(tmp(1,2:3),'_');
end

% get all exp names under the dir, shape is 1*n
all_unique_exp = unique(all_unique_exp);

res = {"exp_name","avg_before","avg_during","avg_after","filter_log"};
% call function and get result
all_matrix=[];

for k = 1 : size(all_unique_exp,2)
    name = all_unique_exp{1,k};
    % get the file name
    for num = 1:4
        fly_file = append('choice_',name,'_N',string(num),'_trackedobjects.csv');
        stimu_file = append('choice_',name,'_stimuli.csv');
        % call function
        try
        [exp_name,avg_before,avg_during,avg_after,filter_res] = ana_pipline(exp_store_dir,fly_file,stimu_file,num,out_dir);
        catch
            exp_name = append(name,'_N',string(num));
            filter_res=deal("No res"); 
            [avg_before,avg_during,avg_after] = deal(nan);
        end

        tmp = {exp_name,avg_before,avg_during,avg_after,filter_res};
        res = [res;tmp];
 
    end
end

% save the result to the current directory
writecell(res,append(out_dir,'summary2.csv'))


function [exp_name,avg_before,avg_during,avg_after,filter_res] = ana_pipline(exp_store_dir,fly_file,stimu_file,num,out_dir);
    %Define ROI
    %clear data_unique;
    
    fly_num = num;

    fly_info_dir = append(exp_store_dir,fly_file);
    stimu_info_dir = append(exp_store_dir,stimu_file);
    name_tmp = strsplit(fly_file,'_');
    exp_name = strjoin(name_tmp(1,2:4),'_');

     %%
    %Import the data 
    data = readmatrix(fly_info_dir);
    %%
    %Import the stimulus
    stimulus=readmatrix(stimu_info_dir);
    %%
    %%%%%%%%%%%%%%%%%%%%%%Load Data Section Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %trimming the data using the ROI
    x_position_raw=data(:,2);
    y_position_raw=data(:,3);
    
    if fly_num==1
    y_position_raw(y_position_raw<=239-4 | y_position_raw>=827-4)=nan;
    x_position_raw(x_position_raw<=310-4 | x_position_raw>=490-4)=nan;
    
    elseif fly_num==2
    y_position_raw(y_position_raw<=243-4 | y_position_raw>=819-4)=nan;
    x_position_raw(x_position_raw<=523-4 | x_position_raw>=700-4)=nan;
    
    elseif fly_num==3
    y_position_raw(y_position_raw<=244-4 | y_position_raw>=820-4)=nan;
    x_position_raw(x_position_raw<=734-4 | x_position_raw>=910-4)=nan;
    
    elseif fly_num==4
    y_position_raw(y_position_raw<=245-4 | y_position_raw>=822-4)=nan;
    x_position_raw(x_position_raw<=944-4 | x_position_raw>=1112-4)=nan;
    end
    
    data(:,2)=x_position_raw;
    data(:,3)=y_position_raw;

    filter_res="";
    if all(data == 0 | isnan(data))
        filter_res=append(filter_res,"empty stuff\n");
        avg_before = nan;
        avg_during = nan;
        avg_after = nan;
        return
    end
    
    
    time=data(:,1);
    x_position=data(:,2);
    y_position=data(:,3);
    original_x_position=x_position;
    original_y_position=y_position;
    %%
    % Calculate speed for each point
    dx = diff(x_position);
    dy = diff(y_position);
    dt = diff(time); % Calculate the time differences

    % Ensure dt is the same size as dx and dy
    dt = dt(1:length(dx));
    
    % Calculate speeds
    speeds = sqrt(dx.^2 + dy.^2) ./ dt; % Speeds at each point

    % Determine threshold for jumping dots
    k = 5; % Constant to be adjusted based on your data
    meanSpeed = nanmean(speeds); % Using nanmean to ignore NaNs
    stdSpeed = nanstd(speeds); % Using nanstd to ignore NaNs
    jumpThreshold = meanSpeed + k * stdSpeed; % Dynamic threshold

    % Detect and handle jumps
    for i = 1:length(speeds)
        if speeds(i) > jumpThreshold
            % Set the next position to NaN to indicate a jump
            if i+1 <= length(x_position)
                x_position(i + 1) = NaN;
                y_position(i + 1) = NaN;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%Filter Section Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Replace the NaN Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the indices of the NaN values in the x_position and y_position vectors
    nan_indices_x = find(isnan(x_position));
    nan_indices_y = find(isnan(y_position));

    % Sort the nan_indices_x and nan_indices_y vectors in ascending order
    nan_indices_x = sort(nan_indices_x);
    nan_indices_y = sort(nan_indices_y);

    % Interpolate NaN values in x_position
    nanIndices = isnan(x_position);
    if any(nanIndices) % Check if there are any NaN values
        x_position(nanIndices) = interp1(find(~nanIndices), x_position(~nanIndices), find(nanIndices), 'linear', 'extrap');
    end

    % Interpolate NaN values in y_position
    nanIndices = isnan(y_position);
    if any(nanIndices) % Check if there are any NaN values
        y_position(nanIndices) = interp1(find(~nanIndices), y_position(~nanIndices), find(nanIndices), 'linear', 'extrap');
    end

    %%
    %replace the data back with x_position
    data_unique(:,2)=x_position;
    data_unique(:,3)=y_position;
    data_unique(:,1)=time;

    %%
    %replace the first jumping point with the a previous valid data

    [~, idx] = unique(data_unique(:,1),'last');
    data_unique = data_unique(idx,:);
    x_position=data_unique(:,2);
    y_position=data_unique(:,3);
    time=data_unique(:,1);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute the velocity
    velocity = (sqrt(diff(x_position).*diff(x_position)+diff(y_position).*diff(y_position))./diff(time))/4;%unit conversion 1mm=4px
    velocity(isinf(velocity))=nan;
    upwind_velocity=(-diff(y_position))./diff(time);%unit conversion 1mm=4px the y axis is inverted
    upwind_velocity=upwind_velocity/4;
    crosswind_velocity=diff(x_position)./diff(time);
    crosswind_velocity=crosswind_velocity/4;

    %%%%%%%%%%%%%%%%%%%%%%%%%%Stimulus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %Matching the stimulus with the time
    stimulus_num=size(stimulus,1);%size(frame) of the stimulus
    time_num=size(time,1);%size(frame) of the total time
    time_stimulus=zeros(size(velocity,1),1);%initialize the time
    for i =1:stimulus_num
       for j=1:time_num
            if(time(j,1)>=stimulus(i,1)&&time(j,1)<=stimulus(i,2))
                time_stimulus(j,1)=1;%time that within the stimulus become 1  
            end
       end
    end
    %Output the time_stimulus with stimulus labeled 1 and non-stimulus with 0
    %%
    %%%%%%%%%%%%%%%%%%%1.Compute the individual speed(scalar)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % find the average of stimulus velocity and non-stimulus velocity
    
    timeframe_stimulus_on=find(time_stimulus==1);
    individual_stimulus_duration=10;
    pre_stimulus_duration=20;
    post_stimulus_duration=40;
    
    time_dot_on_stimuli = find(abs(diff(timeframe_stimulus_on))~=1);
    
    time_dot_on_stimuli = [1;time_dot_on_stimuli+1];
    %time point that stimulus is turning on
    real_time_converter=0.1;

    during_odor_on=timeframe_stimulus_on(time_dot_on_stimuli);

    before_odor=during_odor_on-pre_stimulus_duration/real_time_converter;

    during_odor_off=timeframe_stimulus_on(time_dot_on_stimuli)+individual_stimulus_duration/real_time_converter-1;

    after_odor=during_odor_off+(post_stimulus_duration)/real_time_converter;
    %% 
    %This Part is to use to solve a specific error of tracking, rarely the
    %tracking error is too huge will result a pre-filtered points decreasing,
    %result in a negative time, if that occurs,we discard the data point to
    %ensure the smooth automation of data process

    if any(during_odor_on < 0)
        filter_res=append(filter_res,"discard,tracking error\n");
        avg_before = nan;
        avg_during = nan;
        avg_after = nan;
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%Here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

    %%
    %then let's get the upwind_velocity of before after and during in a array
    %that composes individual trail
    %This Part is to use to solve a specific error of tracking, rarely the
    %tracking error is too huge will result a pre-filtered points decreasing,
    %result in a negative time, if that occurs,we discard the data point to
    %ensure the smooth automation of data process
    %%
    trial_info_matrix=horzcat(before_odor,during_odor_on-1,during_odor_on,during_odor_off,during_odor_off+1,after_odor);%columns contain before,during and after odor information
    %2-3 seconds
    trial_info_matrix_2_3=horzcat(during_odor_on+20,during_odor_on+30);
   
    %%
    %% 
    %create a segment upwind 
    for i=1:length(trial_info_matrix)
        upwind_velocity_seg_before(:,i)=upwind_velocity(trial_info_matrix(i,1):trial_info_matrix(i,2));
        upwind_velocity_seg_during(:,i)=upwind_velocity(trial_info_matrix(i,3):trial_info_matrix(i,4));
        upwind_velocity_seg_during_2_3(:,i)=upwind_velocity(trial_info_matrix_2_3(i,1):trial_info_matrix_2_3(i,2));
        upwind_velocity_seg_after(:,i)=upwind_velocity(trial_info_matrix(i,5):trial_info_matrix(i,6));
    end

    for i=1:length(trial_info_matrix)
        crosswind_velocity_seg_before(:,i)=crosswind_velocity(trial_info_matrix(i,1):trial_info_matrix(i,2));
        crosswind_velocity_seg_during(:,i)=crosswind_velocity(trial_info_matrix(i,3):trial_info_matrix(i,4));
        crosswind_velocity_seg_during_2_3(:,i)=crosswind_velocity(trial_info_matrix_2_3(i,1):trial_info_matrix_2_3(i,2));
        crosswind_velocity_seg_after(:,i)=crosswind_velocity(trial_info_matrix(i,5):trial_info_matrix(i,6));
    end
    %%
    %create the same thing for ground speed
    for i=1:length(trial_info_matrix)
        velocity_seg_before(:,i)=velocity(trial_info_matrix(i,1):trial_info_matrix(i,2));
        velocity_seg_during(:,i)=velocity(trial_info_matrix(i,3):trial_info_matrix(i,4));
        velocity_seg_after(:,i)=velocity(trial_info_matrix(i,5):trial_info_matrix(i,6));
    end
  
    %%
    %create the same thing for x,y coordinate
    for i=1:length(trial_info_matrix)
        x_position_seg_before(:,i)=x_position(trial_info_matrix(i,1):trial_info_matrix(i,2));
        x_position_seg_during(:,i)=x_position(trial_info_matrix(i,3):trial_info_matrix(i,4));
        x_position_seg_after(:,i)=x_position(trial_info_matrix(i,5):trial_info_matrix(i,6));
    end
    for i=1:length(trial_info_matrix)
        y_position_seg_before(:,i)=y_position(trial_info_matrix(i,1):trial_info_matrix(i,2));
        y_position_seg_during(:,i)=y_position(trial_info_matrix(i,3):trial_info_matrix(i,4));
        y_position_seg_after(:,i)=y_position(trial_info_matrix(i,5):trial_info_matrix(i,6));
    end
    
    %%

    for i=1:length(before_odor)
    avg_upwind_velocity_before(i,1)=mean(upwind_velocity_seg_before(:,i));
    avg_upwind_velocity_during(i,1)=nanmean(upwind_velocity_seg_during(:,i));
    avg_upwind_velocity_during_2_3(i,1)=nanmean(upwind_velocity_seg_during_2_3(:,i));
    avg_upwind_velocity_after(i,1)=mean(upwind_velocity_seg_after(:,i));
    end
    %do the same thing for ground speed
    
    for i=1:length(before_odor)
    avg_velocity_before(i,1)=mean(velocity_seg_before(:,i));
    avg_velocity_during(i,1)=mean(velocity_seg_during(:,i));
    avg_velocity_after(i,1)=mean(velocity_seg_after(:,i));
    end

    %%
    %create a filter where pre_odor speed<1mm trials are filtered out
    avg_velocity_total=avg_velocity_before;
    avg_upwind_velocity_after(avg_velocity_total<1)=nan;
    avg_upwind_velocity_before(avg_velocity_total<1)=nan;
    avg_upwind_velocity_during(avg_velocity_total<1)=nan;
    filter_trial_criteria = ones(stimulus_num,1);
    filter_trial_criteria(avg_velocity_total<1)=0;
    %%
    filter_state=0;
    if nanmean(velocity)<1
        filter_res = append(filter_res,"discard,moving speed\n");
        filter_state=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%Here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if nansum(sqrt(diff(x_position).^2+diff(y_position).^2))/4<25
        filter_res = append(filter_res,"discard,moving distance\n");
        filter_state=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%Here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if sum(isnan(avg_upwind_velocity_during))>0.5*length(avg_upwind_velocity_during)
         filter_res = append(filter_res,"discard,can't represent individual moving speed");
         filter_state=1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%Here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %we concatanate things together
    x_position_seg_during(:,filter_trial_criteria == 0)=[];
    y_position_seg_during(:,filter_trial_criteria == 0)=[];
    upwind_trials = vertcat(upwind_velocity_seg_before, upwind_velocity_seg_during, upwind_velocity_seg_after);
    upwind_trials(:, filter_trial_criteria == 0) = [];
    raw_position=horzcat(x_position_seg_during,y_position_seg_during)% so for here the data strcuture would
    %%
    valid_trial_num = size(upwind_trials, 2);
    if filter_state==0
        disp(valid_trial_num);
        average_matrix = sum(upwind_trials,2)/valid_trial_num;
        %writematrix(raw_position,append(out_dir,exp_name,'_pos.csv'));
        writematrix(raw_position,append(out_dir,exp_name,'_pos.csv'))
        avg_before = mean(average_matrix(1:200));
        avg_during = mean(average_matrix(220:230));
        avg_after= mean(average_matrix(301:700));
        %average_matrix(:,2)=valid_trial_num;

    end
    
    %upwind_trials(abs(upwind_trials)>1)=nan;
    %upwind_trials(:,all(isnan(upwind_trials)))=[];
    %writematrix(upwind_trials,append(out_dir,exp_name,'_all.csv'));
    %smooth_upwind_trials=smoothdata(upwind_trials);
    %plot(smooth_upwind_trials)

end
