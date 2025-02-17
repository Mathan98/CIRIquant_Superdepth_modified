import subprocess
import sys
import os


def fastq_cat(input_files, threads, output_file):
    """
    Concatenate multiple FASTQ files into a single output file using Bash commands.
    Handles a mix of gzipped and non-gzipped files. Proceeds only if the output file has a .gz extension.

    Parameters:
        input_files (list): List of input FASTQ file paths.
        output_file (str): Path to the output FASTQ file. Must end with .gz.

    Raises:
        SystemExit: If the output file does not end with .gz or if any subprocess command fails.
    """
    # Ensure the output file has a .gz extension
    if not output_file.endswith('.gz'):
        sys.stderr.write("Error: Output file {} must have a .gz extension.".format(output_file))
        sys.exit(1)

    # Check if all input files are gzipped
    all_gzipped = all(f.endswith('.gz') for f in input_files)
    # Check if no input files are gzipped
    none_gzipped = all(not f.endswith('.gz') for f in input_files)

    try:
        if all_gzipped:
            ## If all files are gzipped, directly concatenate them
            Superdepth = "cat {} > {}".format(" ".join(input_files), output_file)
            subprocess.check_call(Superdepth, shell=True)

        elif none_gzipped:
            ## If no files are gzipped, concatenate and compress the output
            temp_file = 'temp_concatenated.fastq'

            ## Command line
            Superdepth = "cat {} > {}".format(" ".join(input_files), temp_file)
            subprocess.check_call(Superdepth, shell=True)

            gzip = "pigz -p {} -c {} > {}".format(threads, temp_file, output_file)
            subprocess.check_call(gzip, shell=True)
            os.remove(temp_file)

        else:
            ## If there's a mix of gzipped and non-gzipped files
            decompressed_files = []
            decompressed = []
            try:
                for file in input_files:
                    if file.endswith('.gz'):
                        decompressed_file = file[:-3]  # Remove .gz extension

                        ## Command line
                        gunzip = "gunzip -c {} > {}".format(file, decompressed_file)
                        subprocess.check_call(gunzip, shell=True)
                        
                        decompressed.append(decompressed_file)
                        decompressed_files.append(decompressed_file)
                    else:
                        decompressed_files.append(file)

                ## Concatenate all files and compress the output
                temp_file = 'temp_concatenated.fastq'

                ## Command line
                Superdepth = "cat {} > {}".format(" ".join(decompressed_files), temp_file)
                subprocess.check_call(Superdepth, shell=True)
                gzip = "pigz -p {} -c {} > {}".format(threads, temp_file, output_file)
                subprocess.check_call(gzip, shell=True)

                os.remove(temp_file)

            finally:
                ## Clean up decompressed files
                for file in decompressed:
                    if file.endswith('.fastq') and file not in input_files:
                        os.remove(file)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error during concatenation: {}".format(e))
        sys.exit(1)


