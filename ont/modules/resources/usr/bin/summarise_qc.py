#!/usr/bin/env python3

import locale
locale.setlocale(locale.LC_ALL, '')

"""
Mean read length:                 5,615.4
Mean read quality:                   12.9
Median read length:               3,400.0
Median read quality:                 13.6
Number of reads:                744,583.0
Read length N50:                  9,743.0
STDEV read length:                6,015.2
Total bases:              4,181,132,423.0
"""

import sys


key_values = ["Median read length", "Median read quality", "Number of reads", "Read length N50", "Total bases"]

results = {}

for k in key_values: results[k]=[2**63,0]


for fname in sys.argv[1:]:
   with open(fname) as f:
      for line in f:
         data=line.split(":")
         if len(data) !=2: continue
         k=data[0]
         cvt=float if "quality" in k else lambda x: int(x.replace(".0",""))
         if k in results.keys():
             v=cvt(data[1].replace(",","").strip())
             if v<results[k][0] : results[k][0]=v
             if v>results[k][1] : results[k][1]=v

g=open("lr_qc.tex","w")
g.write("""
Table \\ref{tab:lr-qc-summ} summarises the ONT data from the 10 samples successfully sequenced.

\\begin{table}[htb]
\\begin{tabular}{l r r}
Sequencing QC attribute & Minimum & Maximum \\\\\\hline
""")

for k in key_values:
   g.write(f"{k} & {results[k][0]:n} & {results[k][1]:n}\\\\\n")
g.write("""\\hline

\\end{tabular}
\\caption{The minimum and maximum sequencing attribute values  across all samples is shown (after QC).}
\\label{tab:lr-qc-summ}
\\end{table}
""")
g.close()

   
