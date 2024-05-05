#!/usr/bin/env bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                            SCRIPT SUMMARY                             ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## SCRIPT SUMMARY:
## ~~~~~~~~~~~~~~~
# Running execute_run_ABM

 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              PARAMETERS                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Default parameters (correspond to 1 set of Parameter)


parameters_file=""
output_file=""
output_death_file=""



while [[ $# -gt 1   ]]
do
	key="$1"

	case $key in
		-p|--param)
			param="$2"
			shift # past argument
			;;
		-n|--nbr_runs)
			nbr_runs="$2"
			shift # past argument
			;;
		*)
			# unknown option
			;;
	esac
	shift # past argument or value
done


	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              ABM RUN                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
	
#Création fichier log 
LOG=/home/ajaeger/Documents/githubrepos/Stage/log
OUT=/home/ajaeger/Documents/githubrepos/Stage
SCRIPT=/home/ajaeger/Documents/githubrepos/Stage/execute_run_ABM.jl

if [ ! -d $LOG ] ; then mkdir $LOG ; fi
if [ ! -d $LOG/Param_${param} ] ; then mkdir $LOG/Param_${param} ; fi
if [ -f ${OUT}/Param_${param}/parameters_file_${param}.txt ] ; then parameters_file="${OUT}/Param_${param}/parameters_file_${param}.txt" ; fi


for ((ind=1; ind<=nbr_runs; ind++))
do
	if [ ! -d ${OUT}/Param_${param}/output_run_${ind} ]
	then 

	output_file="${OUT}/Param_${param}/output_run_${ind}"
	output_death_file="${OUT}/Param_${param}/output_death_run_${ind}"

	(julia --project=./MyProject $SCRIPT $parameters_file $output_file $output_death_file) >> $LOG//Param_${param}/run_${ind}  2>&1

	fi
done


#Commande pour run le programme
#cd /home/jaeger/sshfs-path_2/Simulation_3/1_Script
#time(bash run_slim.bash -p 1 -c 0 --Ne 25000 --L 10000 --L_HS 100 --n 14 --Start 7000 --Duree_HS 7000 --Duree_tot 14000) 