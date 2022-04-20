#!/usr/bin/python
"""
new and improved MD simulation submition script generator
replaces md_file_gen_cedar.plx perl code

enhancements:
control over the simulation paramaters including:
tempurature, time between frames, time step, min-heat-eq construction
equilibration time
and more...
benchmarking simulations
benchmarking calculations from simulated files
and more.. hopefully. 
"""

#import needed libraries and functions

#import numpy as np
#import pandas as pd
import os
import sys, getopt
import argparse
import configparser


#modules to load for either CPU or GPU SLURM jobs
cpumodules = 'intel/2019.3 openmpi/4.1.1-gnu uofl/amber/18p17_at19p9_intel20193_openmpi312gnu_mkl'
gpumodules = 'intel/2019.3 openmpi/4.1.1-gnu uofl/amber/18p17_at19p9_intel20193_openmpi312gnu_mkl'


#Create the variable parser for all the arguments to go into
parser = argparse.ArgumentParser()

#injest config file if present in current working directory
#use config file flag
parser.add_argument('--config', type=str, required=False, help="(Optional) Provide your config file")

#Simulation Arguments
parser.add_argument('--prmtop', '-p', type=str, required=False, help="(Required) provide your prmtop file") #make this required
parser.add_argument('--prefix', '-i', type=str, required=False, nargs='?', default='replace_me')
parser.add_argument('--segments', '-n', type=int, required=False, nargs='?',help='number of segments to break simulation into')
parser.add_argument('--nstime', '-t', type=int, required=False, nargs='?', help='simulation time per segment')
parser.add_argument('--totaltime', '--tt', type=int, required=False, nargs='?',help='total simulation time per replica')
#SBATCH arguments
parser.add_argument('--walltime', '-w', type=str, required=False, nargs='?')
parser.add_argument('--ntasks', type=str, required=False, nargs='?')
parser.add_argument('--mem', '--ram', type=str, required=False, nargs='?')
parser.add_argument('--partition', type=str, required=False, nargs='?')
parser.add_argument('--mailtype', type=str, required=False, nargs='?',help='all, error, start, cancel or complete, none')
parser.add_argument('--mailuser', type=str, required=False, nargs='?',help='input your email address')
parser.add_argument('--ntwx', type=float, required=False, nargs='?')
#restraint file needed in the .in files?
parser.add_argument('--restraint', type=str, required=False, nargs='?')
parser.add_argument('--constraint', type=str, required=False, nargs='?')

#Creates a variable that stores all the arguments in an array. 
#this new varaible is callable as args.<argument>
args = parser.parse_args()

#Read in Args/flags to fill variables
prmtop     = args.prmtop
prefix     = args.prefix
nsegments  = args.segments
walltime   = args.walltime
ntasks     = args.ntasks
mempercpu  = args.mem            #also known as RAM
partition  = args.partition
mailtype   = args.mailtype
mailuser   = args.mailuser
nstime     = args.nstime         #length of a single chunk
totaltime  = args.totaltime     #length of full simulation (all chunks)
ntwx       = args.ntwx           #default unless overwritten by config or flags
restraint  = args.restraint      #restraint file requested via flag
constraint = args.constraint     #constraint file requested via flag

#Default .in file values unless otherwise specified by a config file
#these values will be output to a default-config.yml file if no config file is provided
#These are PRODUCTION run Default vaules
nsegments  = 1
dt         = 0.002  #default simulation time step (in ps) default is 2fs 
cut        = 10      #default non-bonded cutoff (angstroms)
ntc        = 2       #qmshake restraint. Set to 2 if dt is 2fs or above 
ntt        = 3       #adaptive Termostat to use. Default of 3 is the langevin thermostat (requires gamma_ln)
gamma_ln   = 3.0     #required for ntt. Collision frequency in ps
ntwr       = 5000000 #steps between restart files being written
ntpr       = 50000   #steps between MM energy written
tempi      = 310.0   #initial tempurature
temp0      = 310.0   #reference tempurure (tempurature to maintain the simulation at)
pres0      = 1.0     #reference pressure (in bars) 1 bar = 0.987 atm
imin       = 0       #run minimization flag. 0 off, 1 ON with spe MM calculation
irest      = 1       #restart simulation? 0= generate velocities randomly, 1=restart from .rst7 file
ntx        = 1       #option to read coordinates from the inpcrd file. read formatted with no initial velocity info
ntp        = 1       #constant pressure dynamics, 1=md with isotropic position scaling
ntb        = 2       #periodic boundry conditions. 0=off, 1=const volume, 2=const pressure.
ig         = -1      #seed for the pseudo-random number generator. 
iwrap      = 1       #coordinate 'wrapping'. 1=molecules are re-centered xyz when they pass through a boundry. (aka. autoimaged) 
taup       = 2.0     #pressure relaxation time (in ps) 
ioutfm     = 1       #format out output coordiantes, 1=binary netCDF trajectory
ntr        = 0
ntf        = 2


# If the config flag is enabled and config file provided 
# config parser will read the config file and overwrite variables
configFile = args.config
if configFile is not None:
    print('Using Config file ', configFile)
    configparser = configparser.RawConfigParser()
    configFilePath = os.path.join(os.getcwd(), configFile)
    configparser.read(configFilePath) ##currently the configparger cannot read .in files correctly. 
    #fill (overwrite) arguments with variables from the config file
    #SLURM Variables
    nsegments  = configparser['SLURM']['nsegments']
    walltime   = configparser['SLURM']['walltime']
    ntasks     = configparser['SLURM']['ntasks']
    mempercpu  = configparser['SLURM']['mempercpu']   #also known as RAM
    partition  = configparser['SLURM']['partition']
    mailtype   = configparser['SLURM']['mailtype']
    mailuser   = configparser['SLURM']['mailuser']
    #PMEMD varialbes
    nstime     = configparser['PMEMD']['nstime']
    ntwx       = configparser['PMEMD']['ntwx']
    dt         = configparser['PMEMD']['dt']
    cut        = configparser['PMEMD']['cut']
    ntc        = configparser['PMEMD']['ntc']
    ntt        = configparser['PMEMD']['ntt']
    gamma_ln   = configparser['PMEMD']['gamma_ln']
    ntwr       = configparser['PMEMD']['ntwr']
    ntpr       = configparser['PMEMD']['ntpr']
    tempi      = configparser['PMEMD']['tempi']
    temp0      = configparser['PMEMD']['temp0']
    pres0      = configparser['PMEMD']['pres0']
    imin       = configparser['PMEMD']['imin']
    irest      = configparser['PMEMD']['irest']
    ntx        = configparser['PMEMD']['ntx']
    ntp        = configparser['PMEMD']['ntp']
    ntb        = configparser['PMEMD']['ntb']
    ig         = configparser['PMEMD']['ig']
    iwrap      = configparser['PMEMD']['iwrap']
    taup       = configparser['PMEMD']['taup']
    ioutfm     = configparser['PMEMD']['ioutfm']
    restraint  = configparser['PMEMD']['restraint']
    constraint = configparser['PMEMD']['constraint']
else: #no config file is present (spit out the variables used to "used config file")
    config = configparser.RawConfigParser()
    config['SLURM'] = {'nsegments':nsegments,'walltime':walltime,'ntasks':ntasks,'mempercpu':mempercpu,'partition':partition,'mailtype':mailtype,'mailuser':mailuser}
    config['PMEMD'] = {'nstime':nstime,'ntwx':ntwx,'dt':dt,'cut':cut,'ntc':ntc,'ntt':ntt,'gamma_ln':gamma_ln,'ntwr':ntwr,'ntpr':ntpr,'tempi':tempi,'temp0':temp0,'pres0':pres0,'imin':imin,'irest':irest,'ntx':ntx,'ntp':ntp,'ntb':ntb,'ig':ig,'iwrap':iwrap,'taup':taup,'ioutfm':ioutfm,'restraint':restraint,'constraint':constraint}
    with open('default-run-config.yml', 'w') as configfile:
        config.write(configfile)
        print('No config file provided, Default config written to disk as default-run-config.yml')

##Conditions and manipulation of variables##
#creation of total ntlimit for each segment.
#user can provide
#   --totaltime
#   --nstime
#   --configfile (with nstime) 
if args.totaltime is None and args.nstime is None and configFile is None: 
	parser.error("--totaltime, --nstime args or configfile needed")
if args.totaltime is not None and args.nstime is not None:
	parser.error("either --totaltime or --nstime are required, not both")
if args.totaltime is not None:
	nstlim = int(args.totaltime * 1000 / dt / nsegments) #total time (ns) converted to ps by x1000. divided by time step and number of segment
else:
    args.nstime is not None
    nstlim = int(args.nstime * 1000 / dt)  #time per segment (ns) converted to ps by x1000 and divided by time step 
#summarize nstlim result and where it came from
if args.totaltime is not None:
	 print('--totaltime used for a frame count of ', nstlim)
else:
    args.nstime is not None 
    print('--nstime used for a frame count of ', nstlim)

#now back calculate what the time per chunk is in ns
nschunk = int(nstlim * dt / 1000) #this is used in the MD_in_file to tell the user how many ns each chunk is targeting. 


#this is the .in file that will be used for MD production
md_in_file = f"""
MD_chunk: {nschunk} ns of MD
 &cntrl
  imin = {imin}, irest = {irest}, ntx = {ntx},
  ntb = {ntb}, pres0 = {pres0}, ntp = {ntp},
  ig={ig}, ioutfm={ioutfm},
  iwrap={iwrap},
  taup = {taup},
  cut = {cut}, ntr = {ntr},
  ntc = {ntx}, ntf = {ntf}, ntb={ntb},
  tempi = {tempi}, temp0 = {temp0},
  ntt = {ntt}, gamma_ln = {gamma_ln},
  nstlim = {nstlim}, dt = {dt},
  ntpr = {ntpr}, ntwx = {ntwx}, ntwr = {ntwr},
/
&end
"""



#Variale place holder for min-heat-eq file
minHeatEqSLURM = f"""
#!/bin/bash
#SBATCH --ntasks={ntasks}           # number of MPI processes
#SBATCH --mem-per-cpu={mempercpu}   # memory; default unit is megabytes
#SBATCH --time={walltime}           # time (DD-HH:MM)
#SBATCH --partition={partition}     # 
#SBATCH --mail-type={mailtype}
#SBATCH --mail-user={mailuser}
# mpirun or srun also work

echo "Current working directory is `pwd`"
echo "------------------"

echo "Starting run at: `date`"

module purge
module load {cpumodules}

srun -n {ntasks} pmemd.mpi -O -i {prefix}_01min.in  -o {prefix}_01min.out  -p {prefix}.prmtop -c {prefix}.inpcrd      -r {prefix}_1min.rst7   -ref {prefix}.inpcrd   
srun -n {ntasks} pmemd.mpi -O -i {prefix}_02min.in  -o {prefix}_02min.out  -p {prefix}.prmtop -c {prefix}_1min.rst7   -r {prefix}_2min.rst7   -ref {prefix}_1min.rst7
srun -n {ntasks} pmemd.mpi -O -i {prefix}_03min.in  -o {prefix}_03min.out  -p {prefix}.prmtop -c {prefix}_2min.rst7   -r {prefix}_3min.rst7   -ref {prefix}_2min.rst7
srun -n {ntasks} pmemd.mpi -O -i {prefix}_04min.in  -o {prefix}_04min.out  -p {prefix}.prmtop -c {prefix}_3min.rst7   -r {prefix}_4min.rst7
srun -n {ntasks} pmemd.mpi -O -i {prefix}_05heat.in -o {prefix}_05heat.out -p {prefix}.prmtop -c {prefix}_4min.rst7   -r {prefix}_05heat.rst7 -ref {prefix}_4min.rst7   -x {prefix}_05heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_06heat.in -o {prefix}_06heat.out -p {prefix}.prmtop -c {prefix}_05heat.rst7 -r {prefix}_06heat.rst7 -ref {prefix}_05heat.rst7 -x {prefix}_06heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_07heat.in -o {prefix}_07heat.out -p {prefix}.prmtop -c {prefix}_06heat.rst7 -r {prefix}_07heat.rst7 -ref {prefix}_06heat.rst7 -x {prefix}_07heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_08heat.in -o {prefix}_08heat.out -p {prefix}.prmtop -c {prefix}_07heat.rst7 -r {prefix}_08heat.rst7 -ref {prefix}_07heat.rst7 -x {prefix}_08heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_09heat.in -o {prefix}_09heat.out -p {prefix}.prmtop -c {prefix}_08heat.rst7 -r {prefix}_09heat.rst7 -ref {prefix}_08heat.rst7 -x {prefix}_09heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_10heat.in -o {prefix}_10heat.out -p {prefix}.prmtop -c {prefix}_09heat.rst7 -r {prefix}_10heat.rst7 -ref {prefix}_09heat.rst7 -x {prefix}_10heat.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_11eq.in   -o {prefix}_11eq.out   -p {prefix}.prmtop -c {prefix}_10heat.rst7 -r {prefix}_11eq.rst7   -ref {prefix}_10heat.rst7 -x {prefix}_11eq.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_12eq.in   -o {prefix}_12eq.out   -p {prefix}.prmtop -c {prefix}_11eq.rst7   -r {prefix}_12eq.rst7   -ref {prefix}_11eq.rst7   -x {prefix}_12eq.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_13eq.in   -o {prefix}_13eq.out   -p {prefix}.prmtop -c {prefix}_12eq.rst7   -r {prefix}_13eq.rst7   -ref {prefix}_12eq.rst7   -x {prefix}_13eq.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_14eq.in   -o {prefix}_14eq.out   -p {prefix}.prmtop -c {prefix}_13eq.rst7   -r {prefix}_14eq.rst7   -ref {prefix}_13eq.rst7   -x {prefix}_14eq.nc
srun -n {ntasks} pmemd.mpi -O -i {prefix}_109q.in   -o {prefix}_109q.out   -p {prefix}.prmtop -c {prefix}_14eq.rst7   -r {prefix}_15eq.rst7   -ref {prefix}_14eq.rst7   -x {prefix}_15eq.nc
exit
"""

concatinateSLURM = f"""
#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=0-03:00
#SBATCH --partition=cpu2019,razi-bf,apophis-bf

echo "Current working directory is `pwd`"
echo "------------------"

echo "Starting run at: `date`"

module load {cpumodules}

srun -n 4 cpptraj.mpi -i {prefix}-cpptraj_combine.in > {prefix}-cpptraj_combine.out
exit
"""

#These are Min-heat-eq run values
#min01 min waters, counter ions, freeze: solute
#MINIMIZATION VARIABLES
imin_yes = 1 #minimization IS active
maxcyc = 4000
ncyc = 1000
ntf = 1
iwrap = 1 #YES, create periodic boundry conditions when iwrap = 1
ntwx = 500
ntwe = 500
ntpr = 10
ntr_on = 1 #restraint mask is applied when ntr = 1
ntr_off = 0 #restraint mask is NOT applied when ntr = 0
restraintmask_solvent = "'!(:WAT,Na+,Cl-)'"      #freeze solute
restraintmask_solute_but_not_h = "':WAT,Na+,Cl- & !@H'"   #freeze solvent, but not hydrogens on solute
restraint_wt= 100.0 #force restraint to act on restraint mask
cut = 10.0 #non-bonded calculation cut off set to xx.x angstroms

#HEATING VARIABLES
imin_no = 0 #minimization is not active, regular MD
nstlim_heat = 20000 #20,000 steps of simulation
dt = 0.001 #time step fo 1 femtosecond
ntx = 1 
irest = 0 
ntpr_heat = 250
ntwr_heat_eq = 1000
ntwx_heat = 250
ntt_heat = 3
ig_heat = -1 
gamma_ln_heat = 1.0
ntp_heat = 0
ntc_heat = 1
ntf_heat = 1
ntr_heat = 1
restraint_wt_heat = 25.0

#EQ VARIABLES
nstlim_eq = 10000


#Minimization .in files start here
min01 = f"""
min01: min waters, counter ions, freeze: solute
&cntrl
  imin = {imin_yes}, maxcyc = {maxcyc}, ncyc = {ncyc},
  ntf = {ntf}, iwrap={iwrap},
  ntwx = {ntwx}, ntwe = {ntwe}, ntpr = {ntpr},
  ntr = {ntr_on},
  restraintmask = {restraintmask_solvent},
  restraint_wt= {restraint_wt},
  cut = {cut},
/
&end
"""

min02 = f"""
min02: min solute hydrogen atoms 
&cntrl
  imin = {imin_yes}, maxcyc = {int(maxcyc/2)}, ncyc = {ncyc},
  ntf = {ntf}, iwrap= {iwrap},
  ntwx = {ntwx}, ntwe = {ntwe}, ntpr = {ntpr},
  ntr = {ntr_on},
  restraintmask = {restraintmask_solute_but_not_h}
  restraint_wt= {restraint_wt},
  cut = {cut},
/
&end	
"""

min03 = f"""
min03: min solute 
&cntrl
  imin = {imin_yes}, maxcyc = {int(maxcyc/2)}, ncyc = {ncyc},
  ntf = {ntf}, iwrap={iwrap},
  ntwx = {ntwx}, ntwe = {ntwe}, ntpr = {ntpr},
  ntr = {ntr_on},
  restraintmask = {restraintmask_solvent},
  restraint_wt= {restraint_wt},
  cut = {cut},
/
&end
"""

min04 = f"""
min04: min everything - no restraints
&cntrl
  imin = {imin_yes}, maxcyc = {int(maxcyc-1000)}, ncyc = {ncyc},
  ntf = {ntf}, iwrap= {iwrap},
  ntwx = {ntwx}, ntwe = {ntwe}, ntpr = {ntpr},
  ntr = {ntr_off},
  cut = {cut},
/ 
&end
"""

#HEATING .in files start here
heat05 = f"""
heat05: heat everything with weak restraint on solute 10->60
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx}, irest=0,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=10.0, temp0=60.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

heat06 = f"""
heat06: heat everything with weak restraint on solute 60->110
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx*5}, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=60.0, temp0=110.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

heat07 = f"""
heat07: heat everything with weak restraint on solute 110->160
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx*5}, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=110.0, temp0=160.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

heat08 = f"""
heat08: heat everything with weak restraint on solute 160->260
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx*5}, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250,
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=160.0, temp0=210.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

heat09 = f"""
heat09: heat everything with weak restraint on solute 210->260
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx*5}, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=210.0, temp0=260.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

heat10 = f"""
heat10: heat everything with weak restraint on solute 260->310
&cntrl
  imin= {imin_no},
  nstlim={nstlim_heat}, dt={dt}, ntx={ntx*5}, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=260.0, temp0=310.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt={restraint_wt_heat},
/
&end
"""

#EQUILIBRATION .in files start here
eq11 = f"""
equilibrate11: restraint at 20
&cntrl
  nstlim={nstlim_eq}, dt={dt*2}, ntx=5, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt=20.0,
/ 
&end
"""

eq12 = f"""
equilibrate12: restraint at 15
&cntrl
  nstlim={nstlim_eq}, dt={dt*2}, ntx=5, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt=15.0,
/ 
&end
"""

eq13 = f"""
equilibrate13: restraint at 10
&cntrl
  nstlim={nstlim_eq}, dt={dt*2}, ntx=5, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt=10.0,
/ 
&end
"""

eq14 = f"""
equilibrate14: restraint at 5
&cntrl
  nstlim={nstlim_eq}, dt={dt*2}, ntx=5, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt=5.0,
/ 
&end
"""

eq15 = f"""
equilibrate14: restraint at 1.5
&cntrl
  nstlim={nstlim_eq}, dt={dt*2}, ntx=5, irest=1,
  ntpr=250, ntwr= {ntwr_heat_eq}, ntwx=250,
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = {cut},
  ntr=1, restraintmask = {restraintmask_solvent}, restraint_wt=1.5,
/
&end
EOF
"""
