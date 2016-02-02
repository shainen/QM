from string import Template
from os import getcwd

project=getcwd().split('/')[-1]

runname=getcwd().split('/')[-2]

qsubfile = Template("""
#!/bin/sh
#PBS -j oe
#PBS -l mem=70000mb
#PBS -V
#PBS -m ea
#PBS -M shainen@gmail.com

RUN_NAME=${rname}
SCRATCH_DIR=/data/$$USER/$$RUN_NAME
LOCAL_DIR=/home/shainen/wspin/build

mkdir -p $$SCRATCH_DIR

# Copy them from this directory to data director
cd $$LOCAL_DIR/$$RUN_NAME
cp -r ${prj} $$SCRATCH_DIR/

# Move to the data directory
cd $$SCRATCH_DIR/

# Run the script
time math -script ${prj}/runQM.wl

# Remove the now-useless files
rm -r ${prj} 
""")


with open("../"+runname+".qsub", "w") as f:
    f.write(qsubfile.substitute(rname=runname,prj=project))

