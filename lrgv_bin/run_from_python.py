import os
from subprocess import Popen, PIPE


exe_path = os.path.join(os.path.dirname(__file__), 'lrgvForMOEAFramework.exe')
args = [exe_path, "-m", "std-io", "-b", "AllDecAll", "-c", "ten-year"]
print(f"Launching: {' '.join(args)}")
p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=os.path.dirname(__file__))

# Send input to the executable immediately after starting it
input_line = b"0.5 0.5 0.5 0.5 1.8 0.3 1.9 0.3\n"
print(f"[Python] Sending input to executable: {input_line.decode().strip()}")
p.stdin.write(input_line)
p.stdin.flush()

# Now read output from the executable
output_line = p.stdout.readline()
if output_line:
    print(f"[Executable] {output_line.decode(errors='replace').rstrip()}")
else:
    print("[Python] No output received after sending input.")

print("[Python] Sending final newline to signal end of input.")
p.stdin.write(b"\n")
p.stdin.flush()

# Optionally, read any remaining output
while True:
    line = p.stdout.readline()
    if not line:
        break
    print(f"[Executable] {line.decode(errors='replace').rstrip()}")

print("[Python] Demo complete. Subprocess finished.")