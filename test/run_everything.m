function run_everything()
   send_message("ciao 1");
   send_file_with_telegram("run_everything.m")
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
