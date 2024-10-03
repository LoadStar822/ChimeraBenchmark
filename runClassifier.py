import subprocess
import csv
import os
import sys
import argparse
from datetime import datetime

# Configuration
SOFTWARES = {
    'chimera': {
        'name': 'Chimera',
        'conda_env': 'chimera_env',
        'command': 'chimera -d {db_path} -i {seq_path} -t {threads} -c {threshold} -o {output}'
    },
    'ganon2': {
        'name': 'Ganon2',
        'conda_env': 'ganon2_env',
        'command': 'ganon classify -d {db_path} -i {seq_path} -t {threads} -c {threshold} -o {output}'
    },
    'kraken2': {
        'name': 'Kraken2',
        'conda_env': 'kraken2_env',
        'command': 'kraken2 --db {db_path} --threads {threads} --confidence {threshold} --output {output} {seq_path}'
    },
    'taxor': {
        'name': 'Taxor',
        'conda_env': 'taxor_env',
        'command': 'taxor classify -d {db_path} -i {seq_path} -t {threads} -c {threshold} -o {output}'
    }
}

DEFAULT_OUTPUT_DIR = "./benchmark_results"
DEFAULT_REPORT_FILE = "benchmark_report.csv"

def parse_arguments():
    parser = argparse.ArgumentParser(description="Benchmark multiple bioinformatics tools.")
    parser.add_argument(
        '--software',
        nargs='+',
        choices=SOFTWARES.keys(),
        required=True,
        help='List of software to benchmark. Choices: chimera, ganon2, kraken2, taxor.'
    )
    parser.add_argument(
        '--db-path',
        required=True,
        help='Path to the database.'
    )
    parser.add_argument(
        '--seq-path',
        required=True,
        help='Path to the sequence file (e.g., FASTA file).'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads to use. Default is 4.'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.7,
        help='Classification threshold (0 to 1). Default is 0.7.'
    )
    parser.add_argument(
        '--output-dir',
        default=DEFAULT_OUTPUT_DIR,
        help=f'Output directory. Default is {DEFAULT_OUTPUT_DIR}.'
    )
    return parser.parse_args()

def run_command(cmd, logfile):
    """
    Run the command and save stderr (including time command info) to logfile.
    """
    with open(logfile, 'w') as f:
        process = subprocess.run(cmd, shell=True, stdout=f, stderr=f)
    return process.returncode

def parse_time_log(logfile):
    """
    Parse the time command output log to extract runtime and peak memory.
    """
    runtime = None
    peak_memory = None
    with open(logfile, 'r') as f:
        for line in f:
            if 'Elapsed (wall clock) time' in line:
                # Example format: Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02:34
                try:
                    time_str = line.split(':', 1)[1].strip()
                    parts = time_str.split(':')
                    parts = [float(p) for p in parts]
                    if len(parts) == 3:
                        runtime = parts[0]*3600 + parts[1]*60 + parts[2]
                    elif len(parts) == 2:
                        runtime = parts[0]*60 + parts[1]
                    else:
                        runtime = float(time_str)
                except:
                    runtime = None
            elif 'Maximum resident set size' in line:
                # Unit: KB
                try:
                    mem_kb = float(line.split(':',1)[1].strip())
                    peak_memory = mem_kb / 1024  # Convert to MB
                except:
                    peak_memory = None
    return runtime, peak_memory

def main():
    args = parse_arguments()

    selected_softwares = [SOFTWARES[sw] for sw in args.software]
    db_path = args.db_path
    seq_path = args.seq_path
    threads = args.threads
    threshold = args.threshold
    output_dir = args.output_dir
    results_dir = os.path.join(output_dir, "results")
    report_file = os.path.join(output_dir, DEFAULT_REPORT_FILE)

    # Prepare output directories
    os.makedirs(results_dir, exist_ok=True)

    # Prepare report file
    with open(report_file, 'w', newline='') as csvfile:
        fieldnames = ['Software', 'ReturnCode', 'Runtime_seconds', 'PeakMemory_MB', 'ResultFile']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for software in selected_softwares:
            name = software['name']
            env = software['conda_env']
            cmd_template = software['command']

            print(f"Running {name}...")

            # Construct output file paths
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = os.path.join(results_dir, f"{name}_result_{timestamp}.txt")
            log_file = os.path.join(results_dir, f"{name}_log_{timestamp}.txt")

            # Build actual command
            cmd = cmd_template.format(
                db_path=db_path,
                seq_path=seq_path,
                threads=threads,
                threshold=threshold,
                output=output_file
            )

            # Use GNU time to measure time and memory
            # Ensure to use the full path to GNU time to avoid conflicts
            time_cmd = f"/usr/bin/time -v conda run -n {env} {cmd}"
            full_cmd = time_cmd

            # Run the command
            return_code = run_command(full_cmd, log_file)

            # Parse the log file
            runtime, peak_memory = parse_time_log(log_file)

            if return_code != 0:
                print(f"{name} failed with return code {return_code}. Check log: {log_file}")
            else:
                print(f"{name} completed. Runtime: {runtime} seconds, Peak Memory: {peak_memory} MB")

            # Write to report
            writer.writerow({
                'Software': name,
                'ReturnCode': return_code,
                'Runtime_seconds': runtime if runtime is not None else "N/A",
                'PeakMemory_MB': peak_memory if peak_memory is not None else "N/A",
                'ResultFile': output_file
            })

    print(f"Benchmarking completed. Report saved at {report_file}")

if __name__ == "__main__":
    main()
