#!/usr/bin/env python3

import os
import sys
import argparse
import shlex
import subprocess


def main():
    head_description = '''Bambino wrapper (with hardcoded parameters) that works with RNAIndel.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-b', '--bam', metavar='FILE', required=True,
                        help='input tumor bam file')
    parser.add_argument('-f', '--fasta', metavar='FILE', required=True,
                        help='reference genome FASTA file')
    parser.add_argument('-d', '--dbsnp', metavar='FILE', required=True,
                        help='indels on dbSNP database (https://www.ncbi.nlm.nih.gov/snp) in vcf format')
    parser.add_argument('-o', '--output-file', metavar='FILE', required=True,
                        help='output file with Bambino indel calls')
    args = parser.parse_args()

    # Add Bambino home dir to CLASSPATH
    bambino_home = os.path.dirname(os.path.realpath(__file__))
    try:
        classpath = os.environ['CLASSPATH']
        os.environ['CLASSPATH'] = '{}/*:{}'.format(bambino_home, classpath)
    except KeyError:
        os.environ['CLASSPATH'] = '{}/*'.format(bambino_home)

    # Unpaired Bambino command
    cmd_str = 'java -Xmx6000m Ace2.SAMStreamingSNPFinder -of {} -fasta {} -min-mapq 1 ' \
              '-optional-tags XT!=R -bam {} -tn N -dbsnp-file {} -min-quality 20 ' \
              '-min-flanking-quality 20 -min-alt-allele-count 3 -min-minor-frequency 0 -broad-min-quality 10 ' \
              '-mmf-max-hq-mismatches 6 -mmf-max-hq-mismatches-xt-u 10 -mmf-min-quality 15 -mmf-max-any-mismatches 6 ' \
              '-unique-filter-coverage 2 -no-strand-skew-filter -illumina-q2 1'.format(args.output_file, args.fasta,
                                                                                       args.bam, args.dbsnp)
    stdout, stderr, return_code = run_shell_command(cmd_str)
    if return_code != 0:
        print("Failed while running unpaired Bambino.", file=sys.stderr)
        print(stderr, file=sys.stderr)
        sys.exit(return_code)
    else:
        print('Bambino was successful.', file=sys.stderr)


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

    stdout = proc.stdout.decode('utf-8')
    stderr = proc.stderr.decode('utf-8')

    return_code = proc.returncode

    return stdout, stderr, return_code


if __name__ == '__main__':
    main()
