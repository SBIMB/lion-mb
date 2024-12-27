#!/usr/bin/env python3

import locale
locale.setlocale(locale.LC_ALL, '')

import re
import sys

def get_id(d):
    m=re.search(r">(\w+)",d)
    if m:
       return m.group(1)
    else:
       return d[1:]

seq_dir       = sys.argv[1]
report        = open(sys.argv[2])
gunc_report   = open(sys.argv[3])
completeness  = float(sys.argv[4])
contamination = float(sys.argv[5])
seq_ids       = open(sys.argv[6],"w")

for line in report:
   gline=gunc_report.readline()
   css=float(gline.split()[7].strip())
   data=line.split()
   if len(data)<3: continue
   the_id = data[0]
   compl  = float(data[1])
   contam = float(data[2])
   if compl>=completeness and contam<=contamination and css<0.4:
      f=open(f"{seq_dir}/{the_id}.fa")
      d=f.readline()
      sid=get_id(d)
      seq_ids.write(f"{sid}\n")
seq_ids.close()

