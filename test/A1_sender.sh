#!/bin/bash

# Set your Telegram token and chat_id
telegram_token='6520786123:AAHQxHr0tLDLn62JpzBO0UTvG8TdZ_Dmcos';
chat_id='-827277696;'

# Zip the folder "A1_test"
zip -r A1_test.zip A1_tests

# Function to send the zip file using Telegram API
function send_telegram_message() {
  local file="$1"
  curl -F document=@"$file" "https://api.telegram.org/bot$telegram_token/sendDocument" -F chat_id="$chat_id"
}

# Send the zip file to the chat
send_telegram_message "A1_test.zip"

# Remove the zip file
rm A1_test.zip


