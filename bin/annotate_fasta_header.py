#!/usr/bin/env python3

import re
import locale
locale.setlocale(locale.LC_ALL, '')

#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_kraken_report(kraken_file):
    """
    Parse Kraken2 report file and return a dictionary mapping sequence IDs to taxon IDs.
    
    Args:
        kraken_file (str): Path to Kraken2 report file
        
    Returns:
        dict: Mapping of sequence IDs to taxon IDs
    """
    taxon_mapping = {}
    try:
        with open(kraken_file, 'r') as f:
            for line in f:
                # Split line and extract sequence ID and taxon ID
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    seq_id = fields[1].strip()
                    taxon_id = fields[2].strip()
                    if taxon_id =="2832643": print(fields)
                    taxon_mapping[seq_id] = taxon_id
    except FileNotFoundError:
        sys.stderr.write(f"Error: Kraken report file '{kraken_file}' not found\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error parsing Kraken report: {str(e)}\n")
        sys.exit(1)
    
    return taxon_mapping

def process_fasta(input_fasta, output_fasta, taxon_mapping):
    """
    Process FASTA file and write annotated sequences to output file.
    
    Args:
        input_fasta (str): Path to input FASTA file
        output_fasta (str): Path to output FASTA file
        taxon_mapping (dict): Mapping of sequence IDs to taxon IDs
    """
    try:
        with open(input_fasta, 'r') as in_f, open(output_fasta, 'w') as out_f:
            current_header = ''
            sequence_lines = []
            
            for line in in_f:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if current_header and sequence_lines:
                        write_sequence(out_f, current_header, sequence_lines, taxon_mapping)
                    
                    # Start new sequence
                    current_header = line[1:]  # Remove '>' character
                    sequence_lines = []
                elif line:
                    sequence_lines.append(line)
            
            # Process last sequence
            if current_header and sequence_lines:
                write_sequence(out_f, current_header, sequence_lines, taxon_mapping)
                
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input FASTA file '{input_fasta}' not found\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error processing FASTA file: {str(e)}\n")
        sys.exit(1)

def write_sequence(out_file, header, sequence_lines, taxon_mapping):
    """
    Write a single sequence with annotated header to output file.
    
    Args:
        out_file (file): Output file handle
        header (str): Original sequence header
        sequence_lines (list): Lines containing sequence data
        taxon_mapping (dict): Mapping of sequence IDs to taxon IDs
    """
    # Get sequence ID (first word of header)
    seq_id = header.split()[0]
    
    # Get taxon ID from mapping
    taxon_id = taxon_mapping.get(seq_id, "unknown")
    m=re.search(r"(.*)\(taxid (\d+)\)",taxon_id)
    if m:
        cname=m.group(1)
        taxon_id=m.group(2)
    else:
        cname="unknown"
    
    # Write annotated header and sequence
    out_file.write(f">{seq_id}|kraken:taxid|{taxon_id} {cname}\n")
    out_file.write("\n".join(sequence_lines) + "\n")

def main():
    """Main function to process command line arguments and run the program."""
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python script.py <input_fasta> <kraken_report> <output_fasta>\n")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    kraken_report = sys.argv[2]
    output_fasta   = sys.argv[3]
    
    # Parse Kraken report and get taxon mapping
    taxon_mapping = parse_kraken_report(kraken_report)
    
    # Process FASTA file and write output
    process_fasta(input_fasta, output_fasta, taxon_mapping)

if __name__ == "__main__":
    main()

