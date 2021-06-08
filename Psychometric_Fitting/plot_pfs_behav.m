function plot_pfs_behav(directoryname,figuredirectory)
%plot_pfs_behav(directoryname,figuredirectory)
%
%For each animal in a specified directory, this function
%produces fitted psychometric functions (pfs) within
%individual sessions. This function uses psignifit v4. For additional
%information see: https://github.com/wichmann-lab/psignifit/wiki. Threshold
%and slope values returned are for the scaled fits.
%
%
%Written by MLC 11/28/2016.
%---------------------------------------
warning('off','psignifit:ThresholdPCchanged');
set(0,'DefaultTextInterpreter','none');

[options, plotOptions] = setOptions;

%Get a list of .mat files in the directory
[file_list, file_index] = list_files(directoryname);


%For each file...
for which_file = 1:length(file_index)
    
    %Load file
    filename = file_list(file_index(which_file)).name;
    load([directoryname,'/',filename]);
    
    %Remove fitdata field to start fresh
    if isfield(output,'fitdata')
        output = rmfield(output,'fitdata');
    end
    
    %Set value of dprime that we define as threshold
    options.dprimeThresh = 1;
    
    %Clear plots and handle vectors
    f1 = myplot;
    f2 = myplot;
    handles_f1 = [];
    handles_f2 = [];
    
    
    %For each session...
    for which_session = 1:numel(output)
        
        clear data_to_fit;
        
        %Pull out data from a single session
        data_to_fit = output(which_session).trialmat;
        
        %Continue to next file if there's no data
        if isempty(data_to_fit)
            continue
        end
        
        
        %Fit the data
        [options,results,zFA] = find_threshPC(data_to_fit,options);
        
        %Plot the percent correct values and fit, and save handles
        figure(f1);
        s = subplot(3,4,which_session);
        plotPsych(results,plotOptions);
        handles_f1 = [handles_f1;s]; %#ok<*AGROW>
        
        
        
        %------------------------------------------------------
        %Now transform to dprime space
        %------------------------------------------------------
        figure(f2)
        s2 = subplot(3,4,which_session);
        [x,fitted_yes,fitted_dprime,threshold,slope] = ...
            plotPsych_dprime(results,...
            output(which_session).dprimemat,options,plotOptions,zFA);
        hold on;
        
        %Save plot handles
        handles_f2 = [handles_f2;s2];
        
        
        %Save everything to data structure
        d.results = results;
        d.fit_plot.x = x;
        d.fit_plot.pCorrect = fitted_yes;
        d.fit_plot.dprime = fitted_dprime;
        d.threshold = threshold; %scaled
        d.slope = slope; %scaled
        
        output(which_session).fitdata = d;
        

        
    end
    
    
    %Save figures
    linkaxes(handles_f1);
    linkaxes(handles_f2);
    
    for j = 1:2
        if j == 1
            figType = '_perCorrect';
            f = f1;
        else
            figType = '_dprime';
            f = f2;
        end
        
        fname = [file_list(file_index(which_file)).name(1:end-4),figType];
        suptitle(fname(4:end-4))
        set(f,'PaperPositionMode','auto');
        print(f,'-painters','-depsc', [figuredirectory,fname])
    end
    
    close all
    
    
    %Save file
    fname = file_list(file_index(which_file)).name(1:end-4);
    savename = [directoryname,fname];
    save(savename,'output','-append');
    disp(['Fit and threshold saved successfully to ', savename])
    
end

end





