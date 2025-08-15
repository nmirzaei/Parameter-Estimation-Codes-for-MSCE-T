function coreID = getCoreID()
    if isunix
        [~, cmdout] = system('taskset -cp $$');
        coreID = regexp(cmdout, '\d+', 'match');
    elseif ispc
        [~, cmdout] = system('wmic process where processid="$$" get ProcessorAffinity');
        coreID = regexp(cmdout, '\d+', 'match');
    else
        error('Unsupported platform');
    end
end