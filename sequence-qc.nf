

include {nano_qc_plot as raw_plot } from "./modules/nano-qc.nf" \
    addParams(qc_dir: "${params.qc_dir}/raw",status:"raw")
include {qc_summary; nano_qc_plot as qc_plot } from "./modules/nano-qc.nf" \
    addParams(qc_dir: "${params.qc_dir}/filtered",status:"filtered")
