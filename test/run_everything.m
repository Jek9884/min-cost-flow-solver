function run_everything()
    send_message("==== Start Experiments on A1 ====");
    folder = 'A1_tests';
    recursive_read_folder(folder)
    bashScriptPath = 'A1_sender.sh';
    system(['bash ' bashScriptPath]);
    
    send_message("==== Start Experiments on A2 ====");
    folder = 'A2_tests';
    recursive_read_folder(folder)
    bashScriptPath = 'A2_sender.sh';
    system(['bash ' bashScriptPath]);
    

    disp("End.")
    send_message("End.")
end

function recursive_read_folder(folder)
    fileList = dir(folder);
    for i = 1:length(fileList)
        fileName = fileList(i).name;
        
        if strcmp(fileName, '.') || strcmp(fileName, '..')
            continue;
        end

        pathFile = fullfile(folder, fileName);
        
        if isfile(pathFile)
            str = strtrim(fileName); 
            result = endsWith(str, ".m");
            if result
                launch(pathFile,fileName)
            end

        elseif isfolder(pathFile)
            recursive_read_folder(pathFile)
        end
    end
end

function launch(pathFile,fileName)
    disp(">>> Running experiment:"+fileName)
    send_message(">>> Running experiment:"+fileName)
    run(pathFile);
end

function last_message = send_message(text, chat_id)
    telegram_token = '6520786123:AAHQxHr0tLDLn62JpzBO0UTvG8TdZ_Dmcos';
    if nargin < 2
        chat_id = -827277696;
    end

    try
        url = sprintf('https://api.telegram.org/bot%s/sendMessage?chat_id=%d&text=%s', telegram_token, chat_id, urlencode(text));
        options = weboptions('Timeout', 20);
        output = webread(url, options);

        try
            last_message = output.result.message_id;
        catch
            last_message = [];
        end
    catch
        last_message = [];
    end
end
