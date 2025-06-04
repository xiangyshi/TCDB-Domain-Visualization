#!/usr/bin/env python3

import subprocess
import argparse
import os
import json
import sys

class SaierLabSSH:
    def __init__(self):
        self.script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ssh_execute.sh")
        
        # Ensure the script is executable
        if not os.access(self.script_path, os.X_OK):
            subprocess.run(['chmod', '+x', self.script_path])
    
    def run_commands(self, commands):
        """Run a list of commands on SaierLab server."""
        # Create a temporary modified version of the script with custom commands
        temp_script = "/tmp/temp_ssh_execute.sh"
        
        with open(self.script_path, 'r') as f:
            script_content = f.read()
        
        # Format the commands as a bash array
        commands_str = '    "' + '",\n    "'.join(commands) + '"'
        
        # Replace the COMMANDS section in the script
        modified_content = script_content.replace(
            'COMMANDS=(\n    "echo \'Connected to \\$(hostname)\'"\n    "pwd"\n    "ls -la"\n    "echo \'Current date: \\$(date)\'"\n    # Add more commands here as needed\n)',
            f'COMMANDS=(\n{commands_str}\n)'
        )
        
        with open(temp_script, 'w') as f:
            f.write(modified_content)
        
        # Make the temporary script executable
        subprocess.run(['chmod', '+x', temp_script])
        
        # Run the script
        result = subprocess.run([temp_script], capture_output=True, text=True)
        
        # Clean up
        os.remove(temp_script)
        
        return result.stdout, result.stderr, result.returncode
    
    # Common operations as functions
    def check_system_status(self):
        """Check basic system status on SaierLab."""
        commands = [
            "echo 'Connected to $(hostname)'",
            "uptime",
            "free -h",
            "df -h",
            "who"
        ]
        return self.run_commands(commands)
    
    def list_directory(self, directory="~/"):
        """List contents of a directory on SaierLab."""
        commands = [
            f"echo 'Listing contents of {directory}'",
            f"ls -la {directory}"
        ]
        return self.run_commands(commands)
    
    def run_analysis(self, script_path, input_file=None):
        """Run an analysis script on SaierLab."""
        commands = [
            "cd ~/analysis",
            f"echo 'Running analysis: {script_path}'"
        ]
        
        if input_file:
            commands.append(f"python {script_path} {input_file}")
        else:
            commands.append(f"python {script_path}")
            
        commands.append("echo 'Analysis complete'")
        return self.run_commands(commands)
    
    def backup_data(self, source_dir, backup_name=None):
        """Backup data on SaierLab."""
        if not backup_name:
            backup_name = f"backup_$(date +%Y%m%d_%H%M%S)"
        
        commands = [
            f"echo 'Backing up {source_dir} to {backup_name}'",
            f"mkdir -p ~/backups",
            f"tar -czf ~/backups/{backup_name}.tar.gz {source_dir}",
            f"echo 'Backup complete: ~/backups/{backup_name}.tar.gz'"
        ]
        return self.run_commands(commands)
    
    def custom_commands(self, command_list):
        """Run a custom list of commands."""
        return self.run_commands(command_list)

def main():
    parser = argparse.ArgumentParser(description='Execute commands on SaierLab server')
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Status command
    subparsers.add_parser('status', help='Check system status')
    
    # List directory command
    list_parser = subparsers.add_parser('list', help='List directory contents')
    list_parser.add_argument('directory', nargs='?', default='~/', help='Directory to list')
    
    # Run analysis command
    analysis_parser = subparsers.add_parser('analysis', help='Run analysis script')
    analysis_parser.add_argument('script', help='Analysis script to run')
    analysis_parser.add_argument('input', nargs='?', help='Input file for analysis')
    
    # Backup command
    backup_parser = subparsers.add_parser('backup', help='Backup data')
    backup_parser.add_argument('source', help='Source directory to backup')
    backup_parser.add_argument('--name', help='Backup name')
    
    # Custom command
    custom_parser = subparsers.add_parser('custom', help='Run custom commands')
    custom_parser.add_argument('commands', help='JSON string of commands to run')
    
    args = parser.parse_args()
    ssh = SaierLabSSH()
    
    if args.command == 'status':
        stdout, stderr, code = ssh.check_system_status()
        print(stdout)
        if stderr:
            print(f"Errors: {stderr}", file=sys.stderr)
        
    elif args.command == 'list':
        stdout, stderr, code = ssh.list_directory(args.directory)
        print(stdout)
        if stderr:
            print(f"Errors: {stderr}", file=sys.stderr)
        
    elif args.command == 'analysis':
        stdout, stderr, code = ssh.run_analysis(args.script, args.input)
        print(stdout)
        if stderr:
            print(f"Errors: {stderr}", file=sys.stderr)
        
    elif args.command == 'backup':
        stdout, stderr, code = ssh.backup_data(args.source, args.name)
        print(stdout)
        if stderr:
            print(f"Errors: {stderr}", file=sys.stderr)
        
    elif args.command == 'custom':
        try:
            command_list = json.loads(args.commands)
            if not isinstance(command_list, list):
                raise ValueError("Commands must be a JSON array of strings")
            stdout, stderr, code = ssh.custom_commands(command_list)
            print(stdout)
            if stderr:
                print(f"Errors: {stderr}", file=sys.stderr)
        except json.JSONDecodeError:
            print("Error: Commands must be valid JSON", file=sys.stderr)
            sys.exit(1)
    else:
        parser.print_help()

if __name__ == "__main__":
    main() 