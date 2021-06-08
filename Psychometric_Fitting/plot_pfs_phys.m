function plot_pfs_phys(directoryname,figuredirectory,animalIDs)
%plot_pfs_phys(directoryname,figuredirectory,animalIDs)
%
%This function produces fitted psychometric functions (pfs) for each 
%behavioral session. This function uses psignifit v4. For additional 
%information see https://github.com/wichmann-lab/psignifit/wiki. Threshold
%and slope values returned are for the scaled fits.
%
%
%Input variables:
%directoryname: starting directory that contains the data,
%       organized by test day, then by animal
%
%animalIDs: cellstr array of animal IDs to include in analysis
%
%figuredirectory: directory to store figures
%
%Written by ML Caras Dec 5 2016



%Initialize some parameters
warning('off','psignifit:ThresholdPCchanged');
set(0,'DefaultTextInterpreter','none');
[options, plotOptions] = setOptions;


%Get a list of folders in the directory (each folder == one day of data)
[days,dayIndex]= findRealDirs(directoryname);
days = days(dayIndex);


%For each folder (day)...
for which_day = 1:numel(days)
    
    dayname = [directoryname,days(which_day).name,'/'];
    
    %Get a list of folders in the directory (each folder == one type of data)
    [folders,folderIndex]= findRealDirs(dayname);
    folders = folders(folderIndex);
    
    ind = find_index(folders,'combined');
    
    combined_folder = [dayname,folders(ind).name];
    
    %Get a list of folders in the directory (each folder == one animal)
    [animals,animalIndex] = findRealDirs(combined_folder);
    animals = animals(animalIndex);
    
    %For each animal
    for which_animal = 1:numel(animals)
        
        %Define animal subfolder
        animal_folder = [combined_folder,'/',animals(which_animal).name];
        ID = animals(which_animal).name(end-5:end);
        
        %Make sure animal is one of the ones we want to analyze
        if ~any(ismember(animalIDs,ID))
            continue
        end
        
        
        %Get list of files in the directory
        [files,fileIndex] = listFiles(animal_folder,'*.mat');
        files = files(fileIndex);
        
        %For each file...
        for which_file = 1:numel(files)
            filename = files(which_file).name;
            data_file = [animal_folder,'/',filename];
            
            
            %Only bother fitting data from enagaged sessions. (There is no
            %behavior from disengaged sessions)
            if isempty(strfind(filename,'active'))
                continue
            end
            
            %Load file
            load(data_file)
            
            %Remove fitdata field to start fresh
            if isfield(output,'fitdata')
                output = rmfield(output,'fitdata');
            end
            
            
            %Set value of dprime that we define as threshold
            options.dprimeThresh = 1;
            
            %Clear plots and handle vectors
            f1 = myplot;
            handles_f1 = [];

            clear data_to_fit
            
            %Pull out data for fitting
            data_to_fit = output.trialmat;
            
            %Continue to next file if there's no data
            if isempty(data_to_fit)
                continue
            end
            
            
            %Fit the data
            [options,results,zFA] = find_threshPC(data_to_fit,options);
            
            %Plot the percent correct values and fit, and save handles
            figure(f1);
            s = subplot(2,2,1);
            plotPsych(results,plotOptions);
            handles_f1 = [handles_f1;s]; 
            
            
            
            %------------------------------------------------------
            %Now transform to dprime space
            %------------------------------------------------------
            figure(f1)
            s2 = subplot(2,2,2);
            [x,fitted_yes,fitted_dprime,threshold,slope] = ...
                plotPsych_dprime(results,...
                output.dprimemat,options,plotOptions,zFA);
            hold on;
            
            %Save plot handles and link axes
            handles_f1 = [handles_f1;s2]; %#ok<*AGROW>
            linkaxes(handles_f1,'x');
            
            %Save figure
            fname = [filename(1:end-4),'_psychometric_fits'];
            suptitle([ID,' ',days(which_day).name],f1)
            set(f1,'PaperPositionMode','auto');
            print(f1,'-painter','-depsc', [figuredirectory,fname])
            
            close all
            
            
            %Save everything to data structure
            d.results = results;
            d.fit_plot.x = x;
            d.fit_plot.pCorrect = fitted_yes;
            d.fit_plot.dprime = fitted_dprime;
            d.threshold = threshold; %scaled
            d.slope = slope; %scaled
            
            output.fitdata = d;
            
            
            %Save file
            save(data_file,'output','-append');
            disp(['Fit and threshold saved successfully to ', filename])
           
        end
        
    end
    
end


