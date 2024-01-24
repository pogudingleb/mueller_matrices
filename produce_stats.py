from datetime import datetime
import subprocess
import sys
import numpy as np

## Some helper functions

# Function to run a shell command
def run_shell(args):
    return subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')

# Getting the list of methods
def get_methods():
    full_spec = run_shell(['./matrices'])
    methods = []
    for line in full_spec.split("\n"):
        line = line.strip()
        if line != '' and line[0].isdigit():
            methods.append(line.split(':')[1].strip())
    return methods

# Running a single computation and capturing the runtime
def run_computation(method, num_threads):
    full_output = run_shell(['./matrices', str(method), str(num_threads)])
    rtime = 0
    for line in full_output.split("\n"):
        if line.startswith("Checking performed in"):
            rtime = int(line.split()[-2])
    return rtime

#------------------------------------------

max_threads = 4 if len(sys.argv) == 1 else int(sys.argv[1])
num_runs = 10 if len(sys.argv) <= 2 else int(sys.argv[2])

runtime_data = []
methods = get_methods()
for i, m in enumerate(methods):
    runtime_data.append([])
    for j in range(1, max_threads + 1):
        runtime_data[-1].append([])
        for _ in range(num_runs):
            runtime_data[-1][-1].append(run_computation(i, j))

curtime = datetime.now()
fname = f"timings/timings{curtime.year}_{curtime.month:02}_{curtime.day:02}_{curtime.hour:02}_{curtime.minute:02}.md"
with open(fname, "w") as f:
    f.write(f"Runtimes averaged over {num_runs} runs\n")
    f.write(f"Times reported in microseconds in the form `mean` ± `stdev`\n")

    f.write(f"| Method \ Threads | " + " | ".join(map(str, range(1, max_threads + 1))) + "|\n")
    f.write("|----" * (max_threads + 1) + "|\n")
    
    for m, data in zip(methods, runtime_data):
        time_strings = []
        for series in data:
            avg = sum(series) * 1. / len(series)
            sdev = np.std(series)
            time_strings.append(f"{avg:.1f} ± {sdev:.1f}")
        f.write(f"| {m} |" + " | ".join(time_strings) + "|\n")
