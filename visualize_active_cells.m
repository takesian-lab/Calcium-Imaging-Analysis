function visualize_active_cells(block, plane) 
% Find sound responsive cells within a block and call visualize_cell to
% plot only those cells
%
% Argument(s): 
%   block (struct)
% 
% Returns:
% 
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'

%% Magic numbers & Setup

if nargin > 1
    multiplaneData = true;
    planeName = ['plane' num2str(plane)];
    nPlanes = block.ops.nplanes;
elseif nargin == 1 && isfield(block,'MultiplaneData')
    error('Please choose plane number: visualize_active_cells(block,plane)')
else
    multiplaneData = false;
end

STDlevel = 50;
AUC_F_level = 5;
AUC_S_level = 10;
sort_active = 0;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials

disp('===========PARAMETERS===========')
disp(['STD level = ' num2str(STDlevel)])
disp(['AUC level for df/f = ' num2str(AUC_F_level)])
disp(['AUC level for spks = ' num2str(AUC_S_level)])
disp(['Sort active = ' num2str(sort_active)])
       
%% Find significantly responsive cells
stim_protocol = block.setup.stim_protocol;
baseline_length = block.setup.constant.baseline_length; %seconds
framerate = block.setup.framerate;
nBaselineFrames = round(baseline_length*framerate); %frames

%Accommodate multiplane data
if multiplaneData
    cell_number = block.cell_number.(planeName);
    F7_stim = block.aligned_stim.F7_stim.(planeName);
    spks_stim = block.aligned_stim.spks_stim.(planeName);
else
    cell_number = block.cell_number;
    F7_stim = block.aligned_stim.F7_stim;
    spks_stim = block.aligned_stim.spks_stim;
end

stim_v1 = block.parameters.variable1';
stim_v2 = block.parameters.variable2';

%Identify loco trials to remove
remove = [];
if sort_active == 1
    remove = find(block.active_trials == 1); %active trials
elseif sort_active == 2
    remove = find(block.active_trials == 0); %inactive trials
end

%Separate 0dB 'blank' sound trials
%depending on stim type, 0dB trials are stored in variable1 or variable2
switch stim_protocol
    case {3,2} %{'FM','RF'}
        stim_v0 = stim_v2;

    case {10,9,11} %{'NoiseITI', 'water', 'air'}
        stim_v0 = stim_v1;
        stim_v2 = zeros(size(stim_v1));

    case {5,6} %{'SAM', 'SAMfreq'} %NAN trials instead of zeros           
        stim_v0 = stim_v1;
        stim_v0(isnan(stim_v0)) = 0;

    case {1,12} %Noiseburst, spontaneous
        if length(stim_v1) > 1 || length(stim_v2) > 1
            error('Check stim parameters')
        end
        nTrials = length(block.Sound_Time);
        [stim_v1, stim_v2, stim_v0] = deal(ones(nTrials,1));
        
    otherwise
        error('Stim type is currently not compatible with removing 0dB trials')
end

%Separate blank and stim trials
blankTrials = stim_v0 == 0; %0dB trials
stimTrials = ~blankTrials;

%Get list of all possible stim without blanks prior to removing loco trials
unique_stim_v1 = unique(stim_v1(stim_v0 ~= 0));
unique_stim_v2 = unique(stim_v2(stim_v0 ~= 0));

%If RF, order intensities from highest to lowest
if stim_protocol == 2 %RF
    unique_stim_v2 = flipud(unique_stim_v2);
end
    
%Remove loco trials
blankTrials(remove,:) = 0;
stimTrials(remove,:) = 0;

%Get trial indices
blankTrials = find(blankTrials);
stimTrials = find(stimTrials);

%Only keep stim values for non-blank trials
stim_v1 = stim_v1(stimTrials);
stim_v2 = stim_v2(stimTrials);

%Store stim values in case they change for some cells
store_stim_v1 = stim_v1;
store_stim_v2 = stim_v2;
store_stimTrials = stimTrials;
store_blankTrials = blankTrials;
    
responsiveCells_F = zeros(size(cell_number));
responsiveCells_S = zeros(size(cell_number));

for c = 1:size(cell_number,1)

    %when we remove inf below stim might change so refresh it with original stim list
    stim_v1 = store_stim_v1;
    stim_v2 = store_stim_v2;
    stimTrials = store_stimTrials;
    blankTrials = store_blankTrials;

    F7 = squeeze(F7_stim(c,:,:));
    F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %compute df/f: (total-mean)/mean
    spks = squeeze(spks_stim(c,:,:));

    %Remove trials with infinite values [this was a bug in a small number of blocks]
    [inf_rows,~] = find(isinf(F7_df_f));
    remove_inf = unique(inf_rows);
    if ~isempty(remove_inf)
        stim_v1(stimTrials == remove_inf) = [];
        stim_v2(stimTrials == remove_inf) = [];
        stimTrials(stimTrials == remove_inf) = [];
        blankTrials(blankTrials == remove_inf) = [];
    end

    %Separate stim and blank trials
    F = F7_df_f(stimTrials,:);
    F_blanks = F7_df_f(blankTrials,:);
    S = spks(stimTrials,:);
    S_blanks = spks(blankTrials,:);

    %GET AVERAGED AND SMOOTHED RESPONSES
    %check if each condition is active, then concatenate and keep only active conditions
    analyze_by_stim_condition = 1; %I assume we'll almost always want to do this
    if analyze_by_stim_condition

        F_by_Stim = [];
        S_by_Stim = [];

        for v = 1:length(unique_stim_v1)
            for vv = 1:length(unique_stim_v2)
                stim_rows = intersect(find(stim_v1 == unique_stim_v1(v)), find(stim_v2 == unique_stim_v2(vv)));
                
                if isempty(stim_rows) %Some stim combinations might not exist due to loco
                    continue
                end
                
                F_temp = F(stim_rows,:);
                S_temp = S(stim_rows,:);

                [F_active, ~, ~, ~] = checkIfActive_v2(F_temp, nBaselineFrames, STDlevel, AUC_F_level, 0, '');
                [S_active, ~, ~, ~] = checkIfActive_v2(S_temp, nBaselineFrames, STDlevel, AUC_F_level, 0, '');

                if F_active
                    F_by_Stim = [F_by_Stim; F_temp];
                end

                if S_active
                    S_by_Stim = [S_by_Stim; S_temp];
                end
            end
        end

        %Unless none of the stim conditions end up being significant,
        %use significant conditions only for average
        if ~isempty(F_by_Stim) 
            F = F_by_Stim;
        end

        if ~isempty(S_by_Stim)
            S = S_by_Stim;
        end
    end

    avg_F = nanmean(F,1);
    avg_S = nanmean(S,1);

    %CHECK IF RESPONSIVE
    [responsiveCells_F(c), ~, ~, ~] = checkIfActive_v2(avg_F, nBaselineFrames, STDlevel, AUC_F_level, 0, '');
    [responsiveCells_S(c), ~, ~, ~] = checkIfActive_v2(avg_S, nBaselineFrames, STDlevel, AUC_S_level, 0, '');
    
    %PLOT
    if responsiveCells_F(c)
        figure;
        subplot(1,2,1)
        checkIfActive_v2(avg_F, nBaselineFrames, STDlevel, AUC_F_level, 1, 'df/f');
        
        subplot(1,2,2)
        checkIfActive_v2(avg_S, nBaselineFrames, STDlevel, AUC_S_level, 1, 'spikes');
        
        suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cell_number(c))));
    end
end

disp(['Found ' num2str(sum(responsiveCells_F)) ' responsive F traces and ' num2str(sum(responsiveCells_S)) ' responsive spike traces out of ' num2str(length(cell_number)) ' cells'])

responsiveCells = find(responsiveCells_F); %find(or(responsiveCells_F, responsiveCells_S));
responsiveCellNums = cell_number(responsiveCells);

%% Plot responsive cells only

if isempty(responsiveCellNums)
    return
end

if multiplaneData
    visualize_cell(block, responsiveCellNums, plane)
    visualize_cell(block, responsiveCellNums', plane)
else 
    if stim_protocol == 2 %RF
        visualize_cell_AT_v2(block, responsiveCellNums, [2 4])
    else
        visualize_cell(block, responsiveCellNums)
    end
    visualize_cell(block, responsiveCellNums')
end


end %function end

