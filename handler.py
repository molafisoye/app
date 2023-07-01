#!/usr/bin/env python3

import os
import subprocess
import sys
import boto3
import time
import json

def main():
    start = time.time()
    file_name_prefix = os.environ.get('FILENAME')
    subfolder = os.environ.get('SUBFOLDER')
    target = os.environ.get('TARGET')
    send_status("started", subfolder)
    basename = os.path.basename(file_name_prefix)
    nameroot = basename.split('.')[0]˚
    print(f"starting program with inputs {file_name_prefix}")
    copy_s3_input(file_name_prefix)

    command = ["bash lib/run_call_aberrant_genes.sh", "-i", file_name_prefix, "-o", ".", "-p", nameroot, "-t", target]
    try_run(command)

    subsampled_output = f"{nameroot}_subsampled_{stats}"

    upload_output(f'/app/{subsampled_output}.csv', subfolder)
    upload_output(f'/app/{subsampled_output}.gc_totals.csv', subfolder)
    upload_output(f'/app/normalized.txt', subfolder)

    end = time.time()
    send_status("completed", subfolder, end-start)
    print("The time of execution of secondary pipeline is :", (end-start), "s")

def send_status(status, subfolder, delta=0):
    try:
        status = {
            "jobName":"star",
            "uuid": subfolder,
            "status": status,
            "time": time.time(),
            "delta": delta
        }
        bodyStr = json.dumps(status)
        sqs = boto3.client('sqs')
        response = sqs.send_message(
            QueueUrl="https://sqs.us-east-1.amazonaws.com/287730706223/task-status.fifo",
            MessageGroupId="task-status",
            MessageBody=bodyStr
        )
        print("Response of sending status for STAR...", response)
    except Exception as e:
        print("Exception while send status of STAR...", e)

def try_run(cmd_input):
    if isinstance(cmd_input, list):
        cmd = " ".join(cmd_input)
    else:
        cmd = cmd_input

    with open("error.log","w") as error_log:
        error_log.write(cmd)
        try:
            subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            message = str(e)
            error_log.write(message)
            sys.exit()

def copy_s3_input(input_filename):
    try:
        s3 = boto3.client("s3")
        print(f"downloading {input_filename}")
        s3.download_file(Bucket="gex-fargate-bucket", Key=f"fargate-inputs/{input_filename}", Filename=f"/app/{input_filename}")
    except Exception as e:
        print("exception while copy_s3_input function of secondary pipeline", e)

def upload_output(output_filename, subfolder):
    print(f'uploading - {subfolder}/{output_filename}')
    try:
        s3 = boto3.client("s3")
        s3.upload_file(Bucket="gex-fargate-bucket", Key=f"fargate-outputs/{subfolder}{output_filename}", Filename=output_filename) #output_filename_contains the slash
    except Exception as e:
        print("exception while upload_output function of secondary pipeline", e)

if __name__ == '__main__':
    main()
