#!/bin/bash

# SSH Server Connection and Command Execution Script for SaierLab

# Using your existing SSH config for SaierLab (132.239.144.102)
# No need to specify credentials as they're already in your ~/.ssh/config

# Display connection information
echo "Connecting to SaierLab..."

# Commands to execute on the remote server
# These will be executed in sequence
COMMANDS=(
    "echo 'Connected to \$(hostname)'"
    "pwd"
    "ls -la"
    "echo 'Current date: \$(date)'"
    # Add more commands here as needed
)

# Join commands with semicolons
REMOTE_SCRIPT=$(IFS=';'; echo "${COMMANDS[*]}")

# Connect and execute commands using your existing SSH configuration
ssh SaierLab "$REMOTE_SCRIPT"

echo "SSH session completed."