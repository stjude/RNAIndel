import os
import re
import sys
import shlex
import subprocess
from multiprocessing import Pool

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def callindel(bam, fasta, tmp_dir, heap_memory, region, num_of_processes):

    check_java_version()

    # Add Bambino home dir to CLASSPATH
    bambino_home = os.path.dirname(os.path.realpath(__file__))
    try:
        classpath = os.environ["CLASSPATH"]
        os.environ["CLASSPATH"] = "{}/*:{}".format(bambino_home, classpath)
    except KeyError:
        os.environ["CLASSPATH"] = "{}/*".format(bambino_home)

    # Unpaired Bambino command
    base_cmd_str = (
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

    if region:

        outfile = os.path.join(tmp_dir, "outfile.txt")

        tmp = region.split(":")
        chrom = tmp[0]
        tmp2 = tmp[1].split("-")
        start, stop = tmp2[0], tmp2[1]
        cmd_str = base_cmd_str + " -chr {} -start {} -end {} -of {}".format(
            chrom, start, stop, outfile
        )

        stdout, stderr, return_code = run_shell_command(cmd_str)
        check_caller_return(stdout, stderr, return_code, outfile)
    else:
        if num_of_processes > 1:
            cmds_by_chrom = [
                base_cmd_str
                + " -chr {} -of {}".format(
                    chrom, os.path.join(tmp_dir, "chr{}.txt").format(chrom)
                )
                for chrom in CANONICALS
            ]
            outfiles = [
                os.path.join(tmp_dir, "chr{}.txt").format(chrom) for chrom in CANONICALS
            ]
            pool = Pool(num_of_processes)

            caller_returns = pool.map(run_shell_command, cmds_by_chrom)

            pool.close()
            pool.join()

            for rets, outfile in zip(caller_returns, outfiles):
                check_caller_return(rets[0], rets[1], rets[2], outfile)
        else:

            outfile = os.path.join(tmp_dir, "outfile.txt")

            cmd_str = base_cmd_str + " -of {}".format(outfile)
            stdout, stderr, return_code = run_shell_command(cmd_str)

            check_caller_return(stdout, stderr, return_code, outfile)

    # if return_code != 0 or not os.path.isfile(output_file):
    #    print("Failed while calling indels.", file=sys.stderr)
    #    print(stderr, file=sys.stderr)
    #    sys.exit(return_code)
    # else:
    #    if os.stat(output_file).st_size == 0:
    #        print(
    #            "No variants called. Check if the input reference FASTA file is the same file used for mapping.",
    #            file=sys.stderr,
    #        )
    #        sys.exit(1)
    #    else:
    #        print("indel calling completed successfully.", file=sys.stdout)


def check_java_version():
    pattern = '"(\d+\.{0,1}\d+).*"'

    version = subprocess.check_output(
        ["java", "-version"], stderr=subprocess.STDOUT
    ).decode("utf-8")
    java_ver_lst = re.search(pattern, version).groups()[0].split(".")

    if java_ver_lst:
        if java_ver_lst[0] == "1":
            if int(java_ver_lst[1]) < 8:
                print("Java {} detected.".format(".".join(java_ver_lst)))
                print(
                    "RNAIndel requires Java 8 or higher. Please upgrade Java version."
                )
                sys.exit(1)

        elif int(java_ver_lst[0]) >= 9:
            pass
        else:
            print("Failed to check Java version. Continue anyway...")
    else:
        print("Failed to check Java version. Continue anyway...")


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


def check_caller_return(stdout, stderr, return_code, outfile):

    if return_code != 0:
        print("Failed while calling indels.", file=sys.stderr)
        print(stderr, file=sys.stderr)
        sys.exit(return_code)
    else:
        if not os.path.isfile(outfile):
            print("Failed while calling indels.")
            print(stdout, file=sys.stdout)
            print(stderr, file=sys.stderr)
            sys.exit(1)

        # if os.statdd(outfile).st_size == 0:
        #    print(
        #        "No variants called. Check if the input reference FASTA file is the same file used for mapping.",
        #        file=sys.stderr,
        #    )
        #    sys.exit(1)
        # else:
        #    print("indel calling completed successfully.", file=sys.stdout)


def check_file(file_path, file_name):
    if not os.path.isfile(file_path):
        sys.exit("Error: {} Not Found.".format(file_name))
    return file_path
