#!/usr/bin/env python3
import os
import re
import subprocess

def get_git_files():
    result = subprocess.run(['git', 'ls-files'], capture_output=True, text=True)
    return result.stdout.splitlines()

def filter_file_content(content):
    # Regular expression for email addresses
    email_pattern = r'[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}'
    
    # Split content into lines
    lines = content.split('\n')
    
    # Filter out lines containing email addresses
    filtered_lines = [line for line in lines if not re.search(email_pattern, line)]
    
    # Join lines back together
    return '\n'.join(filtered_lines)

def main():
    # Get list of files
    files = get_git_files()
    
    # Process each file
    for file in files:
        if file.endswith(('.py', '.txt', '.md', '.sh')):
            try:
                # Read file content
                with open(file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Filter content
                filtered_content = filter_file_content(content)
                
                # Write filtered content back to file
                with open(file, 'w', encoding='utf-8') as f:
                    f.write(filtered_content)
                
                print(f"Processed: {file}")
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")

if __name__ == '__main__':
    main() 