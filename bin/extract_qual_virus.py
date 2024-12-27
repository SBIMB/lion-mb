#!/usr/bin/env python3

import pandas as pd
import sys

v_qual_permit= ["High-quality","Medium-quality","Complete"]

df=pd.read_csv(sys.argv[1],sep="\t")

v_qual   = df["checkv_quality"].isin(v_qual_permit)

v_contam = df["contamination"]<=5

v_ok     = df['warnings'].isna()
qual = df[ v_qual & v_contam & v_ok ]
qual["contig_id"].to_csv("quality_viruses.txt",index=False,header=False)


