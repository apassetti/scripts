# Explanation:

- head -n 1 Dataset_S1.txt: Extracts the first row of the file.
- awk -F, '{for (i=1; i<=NF; i++) print length($i)}': Uses awk to print the length of each field in the first row.
- -F, tells awk to use a comma as the field separator.
- for (i=1; i<=NF; i++) print length($i): Iterates through each field and prints its length.
- paste -sd, -: Combines the output into a single comma-separated list of field widths.
- sed 's/,/\t/g': Converts commas to tabs.
- column -ts $'\t' -c $widths: Uses the calculated widths to align each column.
- less -S: Allows horizontal scrolling.

This script should format the output so that each column has a width matching the length of the corresponding field in the first row.
