#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shlex
import subprocess
import tkinter as tk
from tkinter import filedialog

def get_file_path():
    # Create a Tkinter root window (it will not be shown)
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    # Open a file dialog and store the selected file path
    file_path = filedialog.askopenfilename(title="Select a file")
    return file_path

# Get the file path
selected_file = get_file_path()

def main():

    # Get the widths of the columns
    command = f"head -n 1 {shlex.quote(selected_file)} | awk -F, '{{for (i=1; i<=NF; i++) print length($i)}}' | paste -sd, -"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    widths = result.stdout.strip()

    if not widths:
        print("Failed to calculate column widths")
        sys.exit()

    # Execute the script
    command = f"sed 's/,/\\t/g' {shlex.quote(selected_file)} | column -ts $'\\t' -c {widths} | less -S"
    os.system(command)
    print("Done!")

if __name__ == "__main__":
    main()
