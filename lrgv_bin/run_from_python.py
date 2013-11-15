from subprocess import *
p = Popen(["lrgvForMOEAFramework.exe","-m","std-io","-b","AllDecAll","-c","ten-year"], stdin=PIPE, stdout=PIPE)
print p.stdout.readline()
print p.stdout.readline()
print p.stdout.readline()
p.stdin.write("0.5 0.5 0.5 0.5 1.8 0.3 1.9 0.3\n")
print p.stdout.readline()

p.stdin.write("\n")