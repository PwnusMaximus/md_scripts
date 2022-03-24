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

#default Variables
appDescription = "A Python re-build of ryans MD creation script. \
    This script is re-implimented fully in python to create min-heat-eq simulation files as well as md_production files"

#Create the variable parser for all the arguments to go into
parser = argparse.ArgumentParser()

#injest config file if present in current working directory
#use config file flag
parser.add_argument('--config', '-i', type=str, required=False, help="(Optional) Provide your config file")





#Simulation Arguments
parser.add_argument('--prmtop', '-p', type=str, required=True, help="provide your prmtop file")
parser.add_argument('--prefix', '-i', type=str, required=False, nargs='?', default='replace_me')
parser.add_argument('--segments', '-n', type=int, required=False, nargs='?', default=1)
#SBATCH arguments
parser.add_argument('--walltime', '-w', type=str, required=False, nargs='?', default='00-24:00')
parser.add_argument('--ntasks', type=str, required=False, nargs='?', default='8')
parser.add_argument('--mem', '--ram', type=str, required=False, nargs='?', default='2024M')
parser.add_argument('--partition', type=str, required=False, nargs='?', default='cpu2019')
parser.add_argument('--mailtype', type=str, required=False, nargs='?', default='none', help='all, error, start, cancel or complete, none')
parser.add_argument('--mailuser', type=str, required=False, nargs='?', help='input your email address')

#.in md arguments
parser.add_argument('--nstime', '-t', type=int, required=False, nargs='?', default=1)
parser.add_argument('--ntwx', type=float, required=False, nargs='?', default=25000)
parser.add_argument('--dt',type=float, required=False, nargs='?', default=0.002)
parser.add_argument('--cut',type=int, required=False, nargs='?', default=10)
parser.add_argument('--ntc',type=int, required=False, nargs='?', default=2)
parser.add_argument('--ntt',type=int, required=False, nargs='?', default=3)
parser.add_argument('--gamma_ln',type=float, required=False, nargs='?', default=3.0)
parser.add_argument('--ntwr',type=int, required=False, nargs='?', default=5000000)
parser.add_argument('--ntpr',type=int, required=False, nargs='?', default=50000)
parser.add_argument('--tempi',type=float, required=False, nargs='?', default=310.0)
parser.add_argument('--temp0',type=float, required=False, nargs='?', default=310.0)
parser.add_argument('--pres0',type=float, required=False, nargs='?', default=1.0)
parser.add_argument('--imin',type=int, required=False, nargs='?', default=0)
parser.add_argument('--irest',type=int, required=False, nargs='?', default=1)
parser.add_argument('--ntx',type=int, required=False, nargs='?', default=5)
parser.add_argument('--ntp',type=int, required=False, nargs='?', default=1)
parser.add_argument('--ntb',type=int, required=False, nargs='?', default=2)
parser.add_argument('--ig',type=int, required=False, nargs='?', default=-1)
parser.add_argument('--iwrap',type=int, required=False, nargs='?', default=1)
parser.add_argument('--taup',type=float, required=False, nargs='?', default=2.0)
parser.add_argument('--ioutfm',type=int, required=False, nargs='?', default=1)
#restraint file needed in the .in files?
parser.add_argument('--restraint', type=str, required=False, nargs='?')
parser.add_argument('--constraint', type=str, required=False, nargs='?')
#email settings for SLURM

#Creates a variable that stores all the arguments in an array. 
#this new varaible is callable as args.<argument>
args = parser.parse_args()

#Convert Args to variables
prmtop     = args.prmtop
prefix     = args.prefix
nSegments  = args.segments
wallTime   = args.walltime
ntasks     = args.ntasks
memPerCpu  = args.mem          #also known as RAM
partition  = args.partition
mailType   = args.mailtype
mailUser   = args.mailuser
nstime     = args.nstime
ntwx       = args.ntwx
dt         = args.dt
cut        = args.cut
ntc        = args.ntc
ntt        = args.ntt
gamma_ln   = args.gamma_ln
ntwr       = args.ntwr
ntpr       = args.ntpr
tempi      = args.tempi
temp0      = args.temp0
pres0      = args.pres0
imin       = args.imin
irest      = args.irest
ntx        = args.ntx
ntp        = args.ntp
ntb        = args.ntb
ig         = args.ig
iwrap      = args.iwrap
taup       = args.taup
ioutfm     = args.ioutfm
restraint  = args.restraint
constraint = args.constraint

#OS variables



# If the config flag is enabled and config file provided 
# config parser will read the config file and overwrite variables
configFile = args.config
if configFile is not None:
    print('Using Config file ', configFile)
    configparser = configparser.RawConfigParser()
    configFilePath = os.path.join(os.getcwd(), configFile)
    configparser.read(configFilePath) ##currently the configparger cannot read .in files correctly. 
    



#construction of useable variables from the input arguments
nstime = args.nstime * 500000
filePrefix = args.prefix

print(nstime)

print(filePrefix)


#Variale place holder for min-heat-eq file
minHeatEqSLURM = f"""
#!/bin/bash
#SBATCH --ntasks={ntasks}           # number of MPI processes
#SBATCH --mem-per-cpu={memPerCpu}   # memory; default unit is megabytes
#SBATCH --time={wallTime}           # time (DD-HH:MM)
#SBATCH --partition={partition}     # 
#SBATCH --mail-type={mailType}
#SBATCH --mail-user={mailUser}
# mpirun or srun also work

echo "Current working directory is `pwd`"
echo "------------------"

echo "Starting run at: `date`"

module purge
module load intel/2019.3 openmpi/4.1.1-gnu uofl/amber/18p17_at19p9_intel20193_openmpi312gnu_mkl

srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_1min_solv_counterions.in -o {prefix}_1min.out -p {prefix}.prmtop -c {prefix}.inpcrd -r {prefix}_1min.rst7 -ref {prefix}.inpcrd
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_2min_solute_hydrogens.in -o {prefix}_2min.out -p {prefix}.prmtop -c {prefix}_1min.rst7 -ref {prefix}_1min.rst7 -r {prefix}_2min.rst7
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_3min_solute.in -o {prefix}_3min.out -p {prefix}.prmtop -c {prefix}_2min.rst7 -ref {prefix}_2min.rst7 -r {prefix}_3min.rst7
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_4min_all.in -o {prefix}_4min.out -p {prefix}.prmtop -c {prefix}_3min.rst7 -r {prefix}_4min.rst7
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5aheat.in -o {prefix}_5aheat.out -p {prefix}.prmtop -c {prefix}_4min.rst7 -r {prefix}_5aheat.rst7 -ref {prefix}_4min.rst7 -x {prefix}_5aheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5bheat.in -o {prefix}_5bheat.out -p {prefix}.prmtop -c {prefix}_5aheat.rst7 -r {prefix}_5bheat.rst7 -ref {prefix}_5aheat.rst7 -x {prefix}_5bheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5cheat.in -o {prefix}_5cheat.out -p {prefix}.prmtop -c {prefix}_5bheat.rst7 -r {prefix}_5cheat.rst7 -ref {prefix}_5bheat.rst7 -x {prefix}_5cheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5dheat.in -o {prefix}_5dheat.out -p {prefix}.prmtop -c {prefix}_5cheat.rst7 -r {prefix}_5dheat.rst7 -ref {prefix}_5cheat.rst7 -x {prefix}_5dheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5eheat.in -o {prefix}_5eheat.out -p {prefix}.prmtop -c {prefix}_5dheat.rst7 -r {prefix}_5eheat.rst7 -ref {prefix}_5dheat.rst7 -x {prefix}_5eheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_5fheat.in -o {prefix}_5fheat.out -p {prefix}.prmtop -c {prefix}_5eheat.rst7 -r {prefix}_5fheat.rst7 -ref {prefix}_5eheat.rst7 -x {prefix}_5fheat.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_6eq_20.in -o {prefix}_6eq.out -p {prefix}.prmtop -c {prefix}_5fheat.rst7 -r {prefix}_6eq.rst7 -ref {prefix}_5fheat.rst7 -x {prefix}_6eq.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_7eq_15.in -o {prefix}_7eq.out -p {prefix}.prmtop -c {prefix}_6eq.rst7 -r {prefix}_7eq.rst7 -ref {prefix}_6eq.rst7 -x {prefix}_7eq.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_8eq_10.in -o {prefix}_8eq.out -p {prefix}.prmtop -c {prefix}_7eq.rst7 -r {prefix}_8eq.rst7 -ref {prefix}_7eq.rst7 -x {prefix}_8eq.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_9eq_5.in -o {prefix}_9eq.out -p {prefix}.prmtop -c {prefix}_8eq.rst7 -r {prefix}_9eq.rst7 -ref {prefix}_8eq.rst7 -x {prefix}_9eq.nc
srun -n {ntasks} \$AMBERHOME/bin/pmemd -O -i {prefix}_10eq_1.5.in -o {prefix}_10eq.out -p {prefix}.prmtop -c {prefix}_9eq.rst7 -r {prefix}_10eq.rst7 -ref {prefix}_9eq.rst7 -x {prefix}_10eq.nc
exit
"""



#to do immediate, 
    # 

 #prebuilt variables
 
        
"""       
#taken from Ryan's Perl Script
 #!/usr/bin/perl
use warnings;
use strict;
use Math::Trig;
use Getopt::Std;
use List::Util qw(sum);
my %options;
use vars qw(%options);
$options{w} = "1-00:00";
$options{t} = "100";
$options{p} = "replace.prmtop";
$options{i} = "replace";
$options{n} = "1";
getopts ("t:w:p:i:n:hmerc", \%options);
my $timelength = $options{t};
my $steps = $timelength*500000;
my $walltime = $options{w};
my $n = $options{n};
my $prmtop = $options{p};
my $i = 1;
my $j = 1;
my $k = 0;
my $o = 1;
my $q = 2;
my $header = $options{i};

if ($options{h}) {
print "\n  This is a script to generate files for running AMBER md calculations on ARC.\n  The files created are controlled by various flags when executing this script.\n  Once the scipt is competed, you can submit the files by bashing the sub_script.sh\n \n  The -e flag generates min_heat_eq files,\n  the -m flag generates md.sh files and creates cpptraj files to\n  image the restart files between each run of md,\n  the -n flag is used to specify the number of md.sh files to be made (default 1),\n  the -c flag generates a cpptraj file to combine all the created trajectories,\n  the -p flag is used to specify the name of the prmtop file to be used,\n  the -i flag specifies the prefix to be used for the names of the created files,\n  the -t flag sets the amount of time (in ns) for the production in file (default 100ns),\n  the -w flag sets the walltime for the production sh files (DD-HH:MM, default 24h)\n  Please note that ARC has a max gpu walltime of 1 day.\n  The in files should still be modified to change the ntwr, ntwx and ntpr options.\n \n  example usage: perl md_file_gen_arc.plx -e -m -n 5 -c -p wat.prmtop -i water -t 15 -w 0-03:00 \n \n"
}

if ($options{m}) {
&prod;
print "\n Please ensure the created md.in file has the desired ntwr, ntwx and ntpr values\n\n";
}

if ($options{e}) {
&min_heat_eq
}

if ($options{r}) {
&rst
}

if ($options{c}) {
&combine
}

sub min_heat_eq {
	open(my $output_fh, '>', "$header-min_heat_eq.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh <<EOF;
\#!/bin/bash
\#SBATCH --ntasks=8               \# number of MPI processes
\#SBATCH --mem-per-cpu=2024M      \# memory; default unit is megabytes
\#SBATCH --time=0-24:00           \# time (DD-HH:MM)
\#SBATCH --partition=cpu2019
\#SBATCH --mail-type=all
\#SBATCH --mail-user=makay.murray\@uleth.ca
\# mpirun or srun also work

echo "Current working directory is `pwd`"
echo "------------------"

echo "Starting run at: `date`"

module purge
module load intel/2019.3 openmpi/4.1.1-gnu uofl/amber/18p17_at19p9_intel20193_openmpi312gnu_mkl

srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_1min_solv_counterions.in -o $header\_1min.out -p $header.prmtop -c $header.inpcrd -r $header\_1min.rst7 -ref $header.inpcrd
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_2min_solute_hydrogens.in -o $header\_2min.out -p $header.prmtop -c $header\_1min.rst7 -ref $header\_1min.rst7 -r $header\_2min.rst7
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_3min_solute.in -o $header\_3min.out -p $header.prmtop -c $header\_2min.rst7 -ref $header\_2min.rst7 -r $header\_3min.rst7
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_4min_all.in -o $header\_4min.out -p $header.prmtop -c $header\_3min.rst7 -r $header\_4min.rst7
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5aheat.in -o $header\_5aheat.out -p $header.prmtop -c $header\_4min.rst7 -r $header\_5aheat.rst7 -ref $header\_4min.rst7 -x $header\_5aheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5bheat.in -o $header\_5bheat.out -p $header.prmtop -c $header\_5aheat.rst7 -r $header\_5bheat.rst7 -ref $header\_5aheat.rst7 -x $header\_5bheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5cheat.in -o $header\_5cheat.out -p $header.prmtop -c $header\_5bheat.rst7 -r $header\_5cheat.rst7 -ref $header\_5bheat.rst7 -x $header\_5cheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5dheat.in -o $header\_5dheat.out -p $header.prmtop -c $header\_5cheat.rst7 -r $header\_5dheat.rst7 -ref $header\_5cheat.rst7 -x $header\_5dheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5eheat.in -o $header\_5eheat.out -p $header.prmtop -c $header\_5dheat.rst7 -r $header\_5eheat.rst7 -ref $header\_5dheat.rst7 -x $header\_5eheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_5fheat.in -o $header\_5fheat.out -p $header.prmtop -c $header\_5eheat.rst7 -r $header\_5fheat.rst7 -ref $header\_5eheat.rst7 -x $header\_5fheat.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_6eq_20.in -o $header\_6eq.out -p $header.prmtop -c $header\_5fheat.rst7 -r $header\_6eq.rst7 -ref $header\_5fheat.rst7 -x $header\_6eq.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_7eq_15.in -o $header\_7eq.out -p $header.prmtop -c $header\_6eq.rst7 -r $header\_7eq.rst7 -ref $header\_6eq.rst7 -x $header\_7eq.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_8eq_10.in -o $header\_8eq.out -p $header.prmtop -c $header\_7eq.rst7 -r $header\_8eq.rst7 -ref $header\_7eq.rst7 -x $header\_8eq.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_9eq_5.in -o $header\_9eq.out -p $header.prmtop -c $header\_8eq.rst7 -r $header\_9eq.rst7 -ref $header\_8eq.rst7 -x $header\_9eq.nc
srun -n 8 \$AMBERHOME/bin/pmemd -O -i $header\_10eq_1.5.in -o $header\_10eq.out -p $header.prmtop -c $header\_9eq.rst7 -r $header\_10eq.rst7 -ref $header\_9eq.rst7 -x $header\_10eq.nc
exit
EOF

	open(my $output_fh2, '>', "$header\_1min_solv_counterions.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh2 <<EOF;
min01: min waters, counter ions, freeze: solute
\&cntrl
  imin = 1, maxcyc = 4000, ncyc = 1000,
  ntf = 1, iwrap=1,
  ntwx = 500, ntwe = 500, ntpr = 10,
  ntr = 1,
  restraintmask = '!(:WAT,Na+,Cl-)',
  restraint_wt=100.0,
  cut = 10.0,
/
\&end
EOF
	open(my $output_fh3, '>', "$header\_2min_solute_hydrogens.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh3 <<EOF;
min02: min solute hydrogen atoms 
\&cntrl
  imin = 1, maxcyc = 2000, ncyc = 1000,
  ntf = 1, iwrap=1,
  ntwx = 500, ntwe = 500, ntpr = 10,
  ntr = 1,
  restraintmask = ':WAT,Na+,Cl- \& !\@H',
  restraint_wt=100.0,
  cut = 10.0,
/
\&end	
EOF
	open(my $output_fh4, '>', "$header\_3min_solute.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh4 <<EOF;
min03: min solute 
\&cntrl
  imin = 1, maxcyc = 2000, ncyc = 1000,
  ntf = 1, iwrap=1,
  ntwx = 500, ntwe = 500, ntpr = 10,
  ntr = 1,
  restraintmask = '!(:WAT,Na+,Cl-)',
  restraint_wt=100.0,
  cut = 10.0,
/
\&end
EOF
	open(my $output_fh5, '>', "$header\_4min_all.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh5 <<EOF;
min04: min everything - no restraints
\&cntrl
  imin = 1, maxcyc = 3000, ncyc = 1000,
  ntf = 1, iwrap=1,
  ntwx = 500, ntwe = 500, ntpr = 10,
  ntr = 0,
  cut = 10.0,
/ 
\&end
EOF
	open(my $output_fh6, '>', "$header\_5aheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh6 <<EOF;
heat05a: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=1, irest=0,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=10.0, temp0=60.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh7, '>', "$header\_5bheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh7 <<EOF;
heat05b: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=60.0, temp0=110.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh8, '>', "$header\_5cheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh8 <<EOF;
heat05c: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=110.0, temp0=160.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh9, '>', "$header\_5dheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh9 <<EOF;
heat05d: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=160.0, temp0=210.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh10, '>', "$header\_5eheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh10 <<EOF;
heat05e: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=210.0, temp0=260.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh11, '>', "$header\_5fheat.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh11 <<EOF;
heat05f: heat everything with weak restraint on solute
\&cntrl
  imin=0,
  nstlim=20000, dt=0.001, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, ig=-1, gamma_ln = 1.0, tempi=260.0, temp0=310.0,
  ntp=0,
  ntc=1, ntf=1, iwrap=1,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=25.0,
/
\&end
EOF
	open(my $output_fh12, '>', "$header\_6eq_20.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh12 <<EOF;
equilibrate06: restraint at 20
\&cntrl
  nstlim=10000, dt=0.002, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=20.0,
/ 
\&end
EOF
	open(my $output_fh13, '>', "$header\_7eq_15.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh13 <<EOF;
equilibrate07: restraint at 15
\&cntrl
  nstlim=10000, dt=0.002, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=15.0,
/ 
\&end
EOF
	open(my $output_fh14, '>', "$header\_8eq_10.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh14 <<EOF;
equilibrate08: restraint at 10
\&cntrl
  nstlim=10000, dt=0.002, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=10.0,
/ 
\&end
EOF
	open(my $output_fh15, '>', "$header\_9eq_5.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh15 <<EOF;
equilibrate09: restraint at 5
\&cntrl
  nstlim=10000, dt=0.002, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250, 
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=5.0,
/ 
\&end
EOF
	open(my $output_fh16, '>', "$header\_10eq_1.5.in") or die "Couldn't open file $header-md$i.sh :$! \n";
	print $output_fh16 <<EOF;
equilibrate10: restraint at 1.5
\&cntrl
  nstlim=10000, dt=0.002, ntx=5, irest=1,
  ntpr=250, ntwr=1000, ntwx=250,
  ntt=3, temp0=310.0, gamma_ln=1.0, ig=-1,
  ntp=1, taup=2.0,
  ntc=2, ntf=2, iwrap=1, ntb=2,
  cut = 10.0,
  ntr=1, restraintmask = '!(:WAT,Na+,Cl-)', restraint_wt=1.5,
/
\&end
EOF
}

sub prod {
	while ($i <= $n) {
	open(my $output_fh, '>', "$header-md$i.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
	if ($i == 1) {
		print $output_fh <<EOF;
\#!/bin/bash
\#SBATCH --gres=gpu:1
\#SBATCH --mem=1G
\#SBATCH --ntasks=1
\#SBATCH --cpus-per-task=1        # CPU cores per MPI process
\#SBATCH --time=$walltime
\#SBATCH --partition=gpu-v100
\#SBATCH --mail-type=fail
\#SBATCH --mail-user=makay.murray\@uleth.ca
\# mpirun or srun also work

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
nvidia-smi

echo "Current working directory is `pwd`"
echo "------------------"



module load uofl/amber/18_cuda

echo "Starting run at: `date`"
cpptraj -p $prmtop -i $header-md$i-image_rst.in > $header-md$i-image_rst.out
srun -n 1 \$AMBERHOME/bin/pmemd.cuda -O -i $header-md.in -o $header-md$i.out -p $prmtop -c $header\_10eq.ncrst -r $header-md$i.rst7 -x $header-md$i.nc
echo "Finished run at: `date`"
exit
EOF
	open(my $output_fh0, '>', "$header-md$i-image_rst.in") or die "Couldn't open file $header-md$i-image_rst.in :$! \n";
	print $output_fh0 <<EOF;
trajin $header\_10eq.rst7
autoimage
trajout $header\_10eq.ncrst
run
quit
EOF
$i++;
$k++;
next;
	}
	if ($i > 1) {
		print $output_fh <<EOF;
\#!/bin/bash
\#SBATCH --gres=gpu:1
\#SBATCH --mem=1G
\#SBATCH --ntasks=1
\#SBATCH --cpus-per-task=1        # CPU cores per MPI process
\#SBATCH --time=$walltime
\#SBATCH --partition=gpu-v100
\# mpirun or srun also work

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
nvidia-smi

echo "Current working directory is `pwd`"
echo "------------------"


module load uofl/amber/18_cuda

echo "Starting run at: `date`"
cpptraj -p $prmtop -i $header-md$i-image_rst.in > $header-md$i-image_rst.out
srun -n 1 \$AMBERHOME/bin/pmemd.cuda -O -i $header-md.in -o $header-md$i.out -p $prmtop -c $header-md$k.ncrst -r $header-md$i.rst7 -x $header-md$i.nc
echo "Finished run at: `date`"
exit
EOF
	open(my $output_fh0, '>', "$header-md$i-image_rst.in") or die "Couldn't open file $header-md$i-image_rst.in :$! \n";
	print $output_fh0 <<EOF;
trajin $header-md$k.rst7
autoimage
trajout $header-md$k.ncrst
run
quit
EOF
$i++;
$k++;
	}
	}
open(my $output_fh2, '>', "$header-md.in") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh2 <<EOF;
test: $timelength ns of MD
 \&cntrl
  imin = 0, irest = 1, ntx = 5,
  ntb = 2, pres0 = 1.0, ntp = 1,
  ig=-1, ioutfm=1,
  iwrap=1,
  taup = 2.0,
  cut = 10.0, ntr = 0,
  ntc = 2, ntf = 2, ntb=2,
  tempi = 310.0, temp0 = 310.0,
  ntt = 3, gamma_ln = 3.0,
  nstlim = $steps, dt = 0.002,
  ntpr = 50000, ntwx = 25000, ntwr = 5000000,
/
\&end

EOF

if ($options{m} && not (defined($options{e})) && not (defined($options{c})) ) {
open(my $output_fh5, '>', "$header-sub_script.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh5 <<EOF;
\#!/bin/bash

\# first job - no dependencies
jid1=\$(sbatch --parsable $header-md$o.sh)
\# next jobs dependent on the previous jobs

EOF
	while ($q <= $n) {
	print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-md$q.sh) \n";
	$o++;
	$q++;
	}
}

if ($options{m} && (defined($options{e})) && not (defined($options{c})) ) {
open(my $output_fh5, '>', "$header-sub_script.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh5 <<EOF;
\#!/bin/bash

\# first job - no dependencies
jid0=\$(sbatch --parsable $header-min_heat_eq.sh)
\# next jobs dependent on the previous jobs
jid1=\$(sbatch --parsable --dependency=afterany:\$jid0 $header-md$o.sh)

EOF
	while ($q <= $n) {
	print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-md$q.sh) \n";
	$o++;
	$q++;
	}
}

if ($options{e} && not (defined($options{m})) ) {
open(my $output_fh5, '>', "$header-sub_script.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh5 <<EOF;
\#!/bin/bash

\# first job - no dependencies
jid0 = \$(sbatch --parsable $header-min_heat_eq.sh)
EOF
}

if ($options{m} && not (defined($options{e})) && (defined($options{c})) ) {
open(my $output_fh5, '>', "$header-sub_script.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh5 <<EOF;
\#!/bin/bash

\# first job - no dependencies
jid1=\$(sbatch --parsable $header-md$o.sh)
\# next jobs dependent on the previous jobs

EOF
	while ($q <= $n) {
	print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-md$q.sh) \n";
	$o++;
	$q++;
	}
print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-cpptraj_combine.sh) \n";
}

if ($options{m} && (defined($options{e})) && (defined($options{c})) ) {
open(my $output_fh5, '>', "$header-sub_script.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh5 <<EOF;
\#!/bin/bash

\# first job - no dependencies
jid0=\$(sbatch --parsable $header-min_heat_eq.sh)
\# next jobs dependent on the previous jobs
jid1=\$(sbatch --parsable --dependency=afterany:\$jid0 $header-md$o.sh)

EOF
	while ($q <= $n) {
	print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-md$q.sh) \n";
	$o++;
	$q++;
	}
print $output_fh5 "jid$q=\$(sbatch --parsable --dependency=afterany:\$jid$o $header-cpptraj_combine.sh) \n";
}
}

sub combine {
open(my $output_fh3, '>', "$header-cpptraj_combine.sh") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh3 <<EOF;
\#!/bin/bash
\#SBATCH --ntasks=4
\#SBATCH --mem-per-cpu=2024M
\#SBATCH --time=0-03:00
\#SBATCH --partition=cpu2019,razi-bf,apophis-bf
nvidia-smi

echo "Current working directory is `pwd`"
echo "------------------"

echo "Starting run at: `date`"

module load uofl/amber/18p17_at19p9_intel20193_openmpi312gnu_mkl

srun -n 8 \$AMBERHOME/bin/cpptraj -i $header-cpptraj_combine.in > $header-cpptraj_combine.out
exit

EOF

open(my $output_fh4, '>', "$header-cpptraj_combine.in") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh4 <<EOF;
parm $prmtop
EOF
	while ($j <= $n) {
	print $output_fh4 "trajin $header-md$j.nc \n";
	$j++;
	}
print $output_fh4 "autoimage\n";
print $output_fh4 "trajout $header-md-combined.nc netcdf \n";

}

sub rst {
	open(my $output_fh, '>', "rst2pdb_$header.in") or die "Couldn't open file $header-md$i.sh :$! \n";
print $output_fh <<EOF;
module load amber\\16
cpptraj
parm $prmtop
trajin \*md\?.rst7
trajin \*md\?\?.rst7
autoimage origin
trajout $header-md.pdb
run
exit
EOF
}
"""
