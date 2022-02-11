function s = settings(experiment,option)

settings_file = str2func(['settings_' experiment]);
if exist(func2str(settings_file),'file')
    s = settings_file(option);
else
    disp(['Error: No settings file associated with experiment ''' experiment ''' could be found.']);
    s = [];
end
    