import subprocess
import csv
import os
import sys
import argparse
from datetime import datetime
import re
from multitax import NcbiTx

# Import the evaluation scripts as modules
import evaluateClassification
import evaluateAbundance

# Configuration: Base database paths for each software
BASE_DB_PATHS = {
    'chimera': '/mnt/sda/tianqinzhong/project/chimera',
    'ganon2': '/mnt/sdb/mnt_sda/tianqinzhong/project/ganon2',
    'ganon': '/mnt/sdb/mnt_sda/tianqinzhong/project/ganon',
    'kraken2': '/mnt/sdb/mnt_sda/tianqinzhong/project/kraken2',
    'taxor': '/mnt/sdb/mnt_sda/tianqinzhong/project/taxor'
}

# Configuration: Software details
SOFTWARES = {
    'chimera': {
        'name': 'Chimera',
        'conda_env': 'chimera',
        'command': 'chimera classify -i {seq_path} -d {db_path}DB.imcf -t {threads} {threshold_option} -o {output}',
        'threshold_format': '-s {threshold}'
    },
    'ganon2': {
        'name': 'Ganon2',
        'conda_env': 'ganon',
        'command': 'ganon classify -d {db_path} -s {seq_path} -t {threads} {threshold_option} -o {output} --verbose --output-all --output-one',
        'threshold_format': '-c {threshold}'
    },

    'ganon': {
        'name': 'Ganon',
        'conda_env': 'ganon',
        'command': 'ganon classify -d {db_path} -s {seq_path} -t {threads} {threshold_option} -o {output} --verbose --output-all --output-one',
        'threshold_format': '-c {threshold}'
    },

    'kraken2': {
        'name': 'Kraken2',
        'conda_env': 'kraken2',
        'command': 'kraken2 --db {db_path} --threads {threads} {threshold_option} --output {output} {extra} {seq_path} --report {report}',
        'threshold_format': '--confidence {threshold}'
    },

    'taxor': {
        'name': 'Taxor',
        'conda_env': 'taxor',
        'command': 'taxor search --index-file {db_path}.hixf --query-file {seq_path} --output-file {output} {threshold_option} --threads {threads}',
        'threshold_format': '--percentage {threshold}'
    },

    'bracken': {
        'name': 'Bracken',
        'conda_env': 'kraken2',
        'command': 'bracken -d {db_path} -i {report} -o {bracken_output}'
    }
}

DEFAULT_OUTPUT_DIR = "./benchmark_results"
DEFAULT_REPORT_FILE = "benchmark_report.csv"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Benchmark multiple bioinformatics tools.")
    parser.add_argument(
        '-s', '--software',
        nargs='+',
        choices=SOFTWARES.keys(),
        required=True,
        help='List of software to benchmark. Choices: chimera, ganon2, ganon, kraken2, taxor, bracken.'
    )
    parser.add_argument(
        '-d', '--db-subpath',
        required=True,
        help='Database subpath to append to each software\'s base database path.'
    )
    parser.add_argument(
        '-f', '--seq-files',
        nargs='+',
        required=True,
        help='Path(s) to the sequence file(s) (e.g., FASTA files). Supports multiple files.'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help='Number of threads to use. Default is 4.'
    )
    parser.add_argument(
        '-c', '--threshold',
        type=float,
        help='Classification threshold (0 to 1). If not provided, threshold options will not be used in commands.'
    )
    parser.add_argument(
        '-o', '--output-dir',
        default=DEFAULT_OUTPUT_DIR,
        help=f'Output directory. Default is {DEFAULT_OUTPUT_DIR}.'
    )
    # Add optional standard files argument
    parser.add_argument(
        '-sf', '--standard-files',
        nargs='+',
        help='Path(s) to the standard result file(s), corresponding to the sequence files.'
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
    user_time = None
    system_time = None
    with open(logfile, 'r') as f:
        for line in f:
            if 'Elapsed (wall clock) time' in line:
                # Use regex to extract the time string
                match = re.search(r'Elapsed \(wall clock\) time.*:\s*(\d+):(\d+):(\d+)', line)
                if match:
                    hours = int(match.group(1))
                    minutes = int(match.group(2))
                    seconds = int(match.group(3))
                    runtime = hours * 3600 + minutes * 60 + seconds
                else:
                    # Try matching format mm:ss
                    match = re.search(r'Elapsed \(wall clock\) time.*:\s*(\d+):(\d+)', line)
                    if match:
                        minutes = int(match.group(1))
                        seconds = int(match.group(2))
                        runtime = minutes * 60 + seconds
            elif 'User time (seconds)' in line:
                # As a fallback, we can sum user and system time
                try:
                    user_time = float(line.split(':', 1)[1].strip())
                except:
                    user_time = None
            elif 'System time (seconds)' in line:
                try:
                    system_time = float(line.split(':', 1)[1].strip())
                except:
                    system_time = None
            elif 'Maximum resident set size' in line:
                # Unit: KB
                try:
                    mem_kb = float(line.split(':', 1)[1].strip())
                    peak_memory = mem_kb / 1024  # Convert to MB
                except:
                    peak_memory = None
    # If runtime is still None, try to use user_time + system_time
    if runtime is None and user_time is not None and system_time is not None:
        runtime = user_time + system_time
    return runtime, peak_memory


def main():
    args = parse_arguments()

    for seq_file in args.seq_files:
        if not os.path.isfile(seq_file):
            raise FileNotFoundError(f"Sequence file does not exist: {seq_file}")

    # Check if standard files exist (if provided)
    if args.standard_files:
        for std_file in args.standard_files:
            if not os.path.isfile(std_file):
                raise FileNotFoundError(f"Standard file does not exist: {std_file}")

    selected_softwares = [SOFTWARES[sw] for sw in args.software]
    db_subpath = args.db_subpath
    seq_files = args.seq_files
    threads = args.threads
    threshold = args.threshold
    output_dir = args.output_dir
    results_dir = os.path.join(output_dir, "results")

    db_name = os.path.basename(db_subpath.rstrip('/'))  # Get the database name from the subpath
    if args.threshold is not None:
        report_file_name = f"benchmark_report_{db_name}_threshold_{args.threshold}.csv"
    else:
        report_file_name = f"benchmark_report_{db_name}.csv"
    report_file = os.path.join(output_dir, report_file_name)

    # Convert sequence files to absolute paths
    seq_files = [os.path.abspath(f) for f in seq_files]

    # Handle standard files
    if args.standard_files:
        if len(args.standard_files) != len(seq_files):
            print("Error: The number of standard files must match the number of sequence files.")
            sys.exit(1)
        standard_files = [os.path.abspath(f) for f in args.standard_files]
    else:
        standard_files = [None] * len(seq_files)

    print("Starting initialize NcbiDatabase...")
    tax = NcbiTx()

    # Prepare output directories
    os.makedirs(results_dir, exist_ok=True)

    # Prepare report file
    with open(report_file, 'w', newline='') as csvfile:
        fieldnames = ['Software', 'Sequence_File', 'ReturnCode', 'Runtime_seconds', 'PeakMemory_MB', 'ResultFile']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Collect result files for evaluation
        evaluation_tasks = []

        for idx, seq_file in enumerate(seq_files):
            # Ensure the sequence file exists
            if not os.path.isfile(seq_file):
                print(f"Sequence file {seq_file} does not exist. Skipping.")
                continue

            seq_filename = os.path.basename(seq_file)
            seq_name, seq_ext = os.path.splitext(seq_filename)
            # Handle .fasta.gz or similar
            if seq_ext == '.gz':
                seq_name, _ = os.path.splitext(seq_name)
                is_gz = True
            else:
                is_gz = False

            standard_file = standard_files[idx]

            for sw_key in args.software:
                software = SOFTWARES[sw_key]
                name = software['name']
                env = software['conda_env']
                cmd_template = software['command']

                # Get base database path and append subpath
                base_db_path = BASE_DB_PATHS.get(name.lower())
                if not base_db_path:
                    print(f"Base database path for {name} not found. Skipping.")
                    continue
                db_path = os.path.join(base_db_path, db_subpath)

                print(f"Running {name} on {seq_file}...")

                # Construct output file paths
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = os.path.join(results_dir, f"{name}_{seq_name}_{timestamp}.txt")
                log_file = os.path.join(results_dir, f"{name}_{seq_name}_{timestamp}_log.txt")

                # Determine extra parameters
                extra = "--gzip-compressed" if (name.lower() == 'kraken2' and is_gz) else ""

                # For Kraken2, need to specify report file
                if name.lower() == 'kraken2':
                    report_file_kraken = os.path.join(results_dir, f"{name}_{seq_name}_{timestamp}.report")
                else:
                    report_file_kraken = ""

                if args.threshold is not None and 'threshold_format' in software:
                    threshold_option = software['threshold_format'].format(threshold=args.threshold)
                else:
                    threshold_option = ''

                # Build actual command
                cmd = cmd_template.format(
                    db_path=db_path,
                    seq_path=seq_file,
                    threads=threads,
                    threshold_option=threshold_option,
                    output=output_file,
                    extra=extra,
                    report=report_file_kraken,
                    bracken_output=os.path.join(results_dir, f"{seq_name}_{timestamp}.bracken"),
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
                    print(f"{name} on {seq_file} failed with return code {return_code}. Check log: {log_file}")
                else:
                    print(f"{name} on {seq_file} completed. Runtime: {runtime} seconds, Peak Memory: {peak_memory} MB")

                if name.lower() == 'chimera':
                    output_file = output_file + '.tsv'
                elif name.lower() in ['ganon', 'ganon2']:
                    output_file = output_file + '.all'

                # Write to report
                writer.writerow({
                    'Software': name,
                    'Sequence_File': seq_file,
                    'ReturnCode': return_code,
                    'Runtime_seconds': runtime if runtime is not None else "N/A",
                    'PeakMemory_MB': peak_memory if peak_memory is not None else "N/A",
                    'ResultFile': output_file
                })

                # Collect evaluation tasks
                if standard_file:
                    if name.lower() == 'bracken':
                        # For Bracken, only perform abundance evaluation
                        evaluation_tasks.append({
                            'software': 'bracken',
                            'result_file': os.path.join(results_dir, f"bracken_{seq_name}_{timestamp}.bracken"),
                            'standard_file': standard_file,
                            'dataset': seq_name,
                            'database': db_subpath
                        })
                    else:
                        evaluation_tasks.append({
                            'software': sw_key,
                            'result_file': output_file,
                            'standard_file': standard_file,
                            'dataset': seq_name,
                            'database': db_subpath
                        })

                # Handle Bracken separately if selected
                if name.lower() == 'kraken2' and 'bracken' in [sw.lower() for sw in args.software]:
                    # Assuming bracken is selected along with kraken2
                    bracken_sw = SOFTWARES['bracken']
                    bracken_cmd_template = bracken_sw['command']
                    bracken_output = os.path.join(results_dir, f"bracken_{seq_name}_{timestamp}.bracken")
                    bracken_log_file = os.path.join(results_dir, f"bracken_{seq_name}_{timestamp}_log.txt")

                    # Build bracken command
                    bracken_cmd = bracken_cmd_template.format(
                        db_path=db_path,
                        report=report_file_kraken,
                        bracken_output=bracken_output
                    )

                    # Use GNU time to measure time and memory for bracken
                    bracken_time_cmd = f"/usr/bin/time -v conda run -n {bracken_sw['conda_env']} {bracken_cmd}"
                    full_bracken_cmd = bracken_time_cmd

                    # Run bracken
                    bracken_return_code = run_command(full_bracken_cmd, bracken_log_file)

                    # Parse bracken log
                    bracken_runtime, bracken_peak_memory = parse_time_log(bracken_log_file)

                    if bracken_return_code != 0:
                        print(
                            f"Bracken on {seq_file} failed with return code {bracken_return_code}. Check log: {bracken_log_file}")
                    else:
                        print(
                            f"Bracken on {seq_file} completed. Runtime: {bracken_runtime} seconds, Peak Memory: {bracken_peak_memory} MB")

                    if name.lower() == 'chimera':
                        output_file = output_file + '.tsv'
                    elif name.lower() in ['ganon', 'ganon2']:
                        output_file = output_file + '.all'

                    # Write bracken result to report
                    writer.writerow({
                        'Software': 'Bracken',
                        'Sequence_File': seq_file,
                        'ReturnCode': bracken_return_code,
                        'Runtime_seconds': bracken_runtime if bracken_runtime is not None else "N/A",
                        'PeakMemory_MB': bracken_peak_memory if bracken_peak_memory is not None else "N/A",
                        'ResultFile': bracken_output
                    })

                    # Collect evaluation task for Bracken
                    if standard_file:
                        evaluation_tasks.append({
                            'software': 'bracken',
                            'result_file': bracken_output,
                            'standard_file': standard_file,
                            'dataset': seq_name,
                            'database': db_subpath
                        })

    print(f"Benchmarking completed. Report saved at {report_file}")

    # Perform evaluation if standard files are provided
    if evaluation_tasks:
        print("Starting evaluation...")

        for task in evaluation_tasks:
            software = task['software']
            result_file = task['result_file']
            standard_file = task['standard_file']
            dataset = task['dataset']
            database = task['database']
            output_prefix = os.path.join(output_dir, f"{software}_{dataset}_{database}")

            if software == 'bracken':
                # Only perform abundance evaluation
                print(f"Evaluating abundance for {software} on dataset {dataset}...")
                evaluateAbundance.evaluate_abundance(
                    standard_file=standard_file,
                    result_files=[result_file],
                    software_list=[software],
                    output_file=f"{output_prefix}_abundance.csv",
                    dataset=dataset,
                    database=database,
                    tax=tax
                )
            else:
                # Perform both classification and abundance evaluation
                print(f"Evaluating classification and abundance for {software} on dataset {dataset}...")
                # Classification evaluation
                evaluateClassification.evaluate_classification(
                    standard_file=standard_file,
                    result_file=result_file,
                    software=software,
                    output_file=f"{output_prefix}_classification.csv",
                    dataset=dataset,
                    database=database,
                    tax=tax
                )
                # Abundance evaluation
                evaluateAbundance.evaluate_abundance(
                    standard_file=standard_file,
                    result_files=[result_file],
                    software_list=[software],
                    output_file=f"{output_prefix}_abundance.csv",
                    dataset=dataset,
                    database=database,
                    tax=tax
                )

        print("Evaluation completed.")


if __name__ == "__main__":
    main()
