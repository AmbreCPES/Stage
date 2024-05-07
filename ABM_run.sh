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
data_overwrite="false"



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
		-steps|--n_steps)
			n_steps="$2"
			shift # past argument
			;;
		-ow|--overwrite)
			data_overwrite="$2"
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
if [ -f ${OUT}/Param_${param}/parameters_file_${param}.txt ]
then 
	parameters_file="${OUT}/Param_${param}/parameters_file_${param}.txt"
else 
	echo -e "$\n{OUT}/Param_${param}/parameters_file_${param}.txt does not exist"
	exit 0
fi


for ((ind=1; ind<=nbr_runs; ind++))
do
	if [ $data_overwrite = "true" ] || [ ! -f ${OUT}/Param_${param}/output_run_${ind} ]
	then 
		
		output_file="${OUT}/Param_${param}/output_run_${ind}"
		output_death_file="${OUT}/Param_${param}/output_death_run_${ind}"

		echo -e "Dead cells information stored in: $output_death_file \nCells information stored in: $output_file \n" >> $LOG//Param_${param}/run_${ind}  2>&1
		less $parameters_file >> $LOG//Param_${param}/run_${ind}  2>&1

		time(julia --project=./MyProject $SCRIPT $n_steps $parameters_file $output_file $output_death_file) >> $LOG//Param_${param}/run_${ind}  2>&1

	else

		echo -e "\n${OUT}/Param_${param}/output_run_${ind} already exists and overwrite not allowed, to overwrite add the option -ow "true" " >> $LOG//Param_${param}/run_${ind}  2>&1
		exit 0
	fi
done


#Commande pour run le programme
#cd /home/jaeger/sshfs-path_2/Simulation_3/1_Script
#bash ABM_run.sh -p 0 -n 1 -steps 100  