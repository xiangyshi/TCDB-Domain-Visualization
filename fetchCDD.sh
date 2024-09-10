#!/bin/bash

# Variables
REMOTE_HOST="SaierLab"
REMOTE_FILE_PATH="/Users/leo/test/e4.cdd"  # Adjust this path if needed
LOCAL_DESTINATION_PATH="/Users/leo/Desktop/saierlab/"

# Step 1: Securely copy the file from the remote machine to the local machine
scp ${REMOTE_HOST}:${REMOTE_FILE_PATH} ${LOCAL_DESTINATION_PATH}

# Optional: Notify the user of completion
echo "File successfully copied to ${LOCAL_DESTINATION_PATH}"
