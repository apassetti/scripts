# scripts
List of useful bioinfo scripts

## gunzip -c data/gencode.v46.annotation.gff3.gz | grep -v "#" | awk '$3 == "gene" && $9 ~ "FKBP1A"' | tee FKBP1A.txt | less -S
This command will:

1. Unzip the GFF3 file.
2. Filter out comment lines (lines starting with `#`).
3. Select only lines where the third column is `"gene"` and the ninth column matches `"FKBP1A"`.
4. Write the output to `FKBP1A.txt` using `tee`.
5. Display the output in the `less` viewer with horizontal scrolling (`-S` option).

You'll see the output in the terminal while it's also saved to `FKBP1A.txt`.

## gunzip -c hg002.wf_snp.vcf.gz| grep -v "#" |  awk '$8 =="P"' | less -S
This command will:

1. Unzip the `hg002.wf_snp.vcf.gz` file.
2. Filter out comment lines (those starting with `#`).
3. Select only lines where the 8th column equals `"P"`.
4. Save the output to `output.txt` using `tee`.
5. Display the output in `less` with horizontal scrolling enabled (`-S` option).

This way, you can view the output on the terminal with `less` while simultaneously saving it to a file.
