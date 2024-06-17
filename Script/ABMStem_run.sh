#!/usr/bin/env bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                            SCRIPT SUMMARY                             ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## SCRIPT SUMMARY:
## ~~~~~~~~~~~~~~~
# Running ABMStem_Run.jl

 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              PARAMETERS                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Default parameters (correspond to 1 set of Parameter)


parameters_file=""
output_file=""
data_overwrite="false"
start_ind=1
end_ind=10


while [[ $# -gt 1   ]]
do
	key="$1"

	case $key in
		-p|--param)
			param="$2"
			shift # past argument
			;;
		-s|--start)
			start_ind="$2"
			shift # past argument
			;;
		-e|--end)
			end_ind="$2"
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
SCRIPT=/home/ajaeger/Documents/githubrepos/Stage/Script/ABMStem_Run.jl

if [ ! -d $LOG ] ; then mkdir $LOG ; fi
if [ ! -d $LOG/Stem_Param_${param} ] ; then mkdir $LOG/Stem_Param_${param} ; fi
if [ -f ${OUT}/Stem_Param_${param}/stem_parameters_file_${param}.txt ]
then 
	parameters_file="${OUT}/Stem_Param_${param}/stem_parameters_file_${param}.txt"
else 
	echo -e "${OUT}/Stem_Param_${param}/stem_parameters_file_${param}.txt"
	exit 0
fi


for ((ind=start_ind; ind<=end_ind; ind++))
do
	if [ $data_overwrite = "true" ] || [ ! -f ${OUT}/Stem_Param_${param}/output_run_${ind} ]
	then 
		test_mean="${ind}"
		histo_file="${OUT}/Stem_Param_${param}/Analysis/histo_act_param_${param}_${ind}"
		output_file="${OUT}/Stem_Param_${param}/Data/output_act_run_${ind}"
        count_file="${OUT}/Stem_Param_${param}/Data/count_act_run_${ind}"
        plt_file="${OUT}/Stem_Param_${param}/Analysis/plt_run_${ind}"
		
		echo -e "\n############ BEGIN LOG ##########\n" >> $LOG/Stem_Param_${param}/run_${ind} 2>&1
		echo -e "Cells information stored in: $output_file \n" >> $LOG/Stem_Param_${param}/run_${ind}  2>&1

		echo -e "\n############ PARAMETER FILE ##########\n" >> $LOG/Stem_Param_${param}/run_${ind} 2>&1
		less $parameters_file >> $LOG/Stem_Param_${param}/run_${ind}  2>&1
		echo -e "\n############ END PARAMETER FILE ##########\n" >> $LOG/Stem_Param_${param}/run_${ind} 2>&1

		time(julia --project=./MyProject $SCRIPT $n_steps $parameters_file $output_file $count_file $histo_file $plt_file $test_mean) >> $LOG/Stem_Param_${param}/run_${ind}  2>&1

	else

		echo -e "\n${OUT}/Stem_Param_${param}/Data/output_run_${ind} already exists and overwrite not allowed, to overwrite add the option -ow "true" " >> $LOG/Stem_Param_${param}/run_${ind}  2>&1
		exit 0
	fi
done


#Commande pour run le programme
#bash ABMStem_run.sh -p 0 -s 1 -e 10 -steps 100 &