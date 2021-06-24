import os
import sys
import shlex
import subprocess
from functools import partial


def callindel(bam, fasta, output_file, heap_memory, region):
    
    # Add Bambino home dir to CLASSPATH
    bambino_home = os.path.dirname(os.path.realpath(__file__))
    try:
        classpath = os.environ["CLASSPATH"]
        os.environ["CLASSPATH"] = "{}/*:{}".format(bambino_home, classpath)
    except KeyError:
        os.environ["CLASSPATH"] = "{}/*".format(bambino_home)

    # Unpaired Bambino command
    cmd_str = (
        "java -Xmx{} Ace2.SAMStreamingSNPFinder -fasta {} -min-mapq 1 "
        "-optional-tags XT!=R -bam {} -tn T -min-quality 20 "
        "-min-flanking-quality 20 -min-alt-allele-count 3 "
        "-min-minor-frequency 0 -broad-min-quality 10 "
        "-mmf-max-hq-mismatches 8 -mmf-max-hq-mismatches-xt-u 10 "
        "-mmf-min-hq-quality 15 -mmf-max-lq-mismatches 8 "
        "-unique-filter-coverage 2 -no-strand-skew-filter -illumina-q2 1 "
        "-poly-x-min-run-length 10 -merge-equivalent-indels -autotune -query-mode".format(
            heap_memory, fasta, bam
        )
    )
    
    cmd_str = cmd_str + " -of {}".format(output_file)

    if region:
        # regex check: expected pattern (chr)[0-9XY]+:[0-9]+-[0-9]+

        tmp = region.split(":") 
        chrom = tmp[0]
        tmp2 = tmp[1].split("-")
        start, stop = tmp2[0], tmp2[1]
        #chrom, start, stop = region[0], region[1], region[2]
        cmd_str = cmd_str + " -chr {} -start {} -end {}".format(chrom, start, stop)

    stdout, stderr, return_code = run_shell_command(cmd_str)

    if return_code != 0 or not os.path.isfile(output_file):
        print("Failed while calling indels.", file=sys.stderr)
        print(stderr, file=sys.stderr)
        sys.exit(return_code)
    else:
        if os.stat(output_file).st_size == 0:
            print(
                "No variants called. Check if the input reference FASTA file is the same file used for mapping.",
                file=sys.stderr,
            )
            sys.exit(1)
        else:
            print("indel calling completed successfully.", file=sys.stdout)


def run_shell_command(command_string):
    """ Executes a command and returns stdout, stderr, return_code.
        Input:
            - command_string: Command to be executed
        Output:
            - stdout: stdout of command as a single string.
            - stderr: stderr of command as a single string.
            - return_code: integer return code of command.
    """
    command = shlex.split(command_string)
    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout = proc.stdout.decode("utf-8")
    stderr = proc.stderr.decode("utf-8")

    return_code = proc.returncode

    return stdout, stderr, return_code


def check_file(file_path, file_name):
    if not os.path.isfile(file_path):
        sys.exit("Error: {} Not Found.".format(file_name))
    return file_path

