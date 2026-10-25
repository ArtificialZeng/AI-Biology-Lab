#!/bin/bash
set -e

# usage
#taskset -c 50-112 bash run_ESM2ne_viral.sh 
########################################################################################################
#                                      Viral sequences
#                                   Debora marks datasets 
########################################################################################################
viral=(
    # Singles
    SARS2_RBD_binding_Starr       SARS2_RBD_expression_Starr    SARS2_BA1_SPIKE_Dadonaite        SARS2_DELTA_SPIKE_Dadonaite   
    SARS2_PRD0038_RBD_Starr       SARS2_PLPRO_activity_Wu       SARS2_PLPRO_abundance_Wu         SARS2_MRPO_Flynn
    SARS2_XBB15_RBD_Taylor                 
    
    IAV_NA_Jiang                  IAV_H1_HA_Wu                  IAV_H1_NP_Doud                   IAV_H5_HA_Dadonaite     
    IAV_H1_HA_Doud                IAV_PA_Wu                     IAV_PB2_Soh                      IAV_RDRP_Li                 
    IAV_H3_NP_Doud                IAV_H3_HA_Lee
    
    CVB3_2A_Alvarez               CVB3_2B_Alvarez               CVB3_2C_Alvarez                  CVB3_3A_Alvarez    
    CVB3_3B_Alvarez               CVB3_3C_Alvarez               CVB3_3D_Alvarez                  CVB3_VP1_Alvarez
    CVB3_VP3_Alvarez              CVB3_POLG_Mattenberger
    
    HIV1_BF520_ENV_Haddox         HIV1_BG505_ENV_Haddox         HIV1_HV1B9_ENV_DuenasDecamp      NIPAH_F_Larsen
    RmYN02_RBD_Starr              RsYN04_RBD_Starr              DENV_POLG_Suphatrakul            LAMBDA_HCP_Tsuboyama 
    LASSA_GP_Carr                 EV_CAPSD_Bakhache             EV_REP_Bakhache                  BPP22_COAT_Tsuboyama          
                                 
    ## Doubles and multiple      
    # PESV_POLG_Tsuboyama     BP434_RPC1_Tsuboyama
    ## Multiple
    #AAV2_CAPSD_Sinai
    )


######################### Define variables ########################
# function to generate a seed
make_seed() {
    local dataset_id=$1
    local rep=$2
    local s="${dataset_id}_${rep}"
    local h=$(echo -n "$s" | md5sum | cut -c1-8)
    printf "%d\n" "0x$h"
}


checkpoint="esm2_t33_650M_UR50D"
strategies=('pool_split' 'site_split')
source="viral"

for rep in 1 2 3
do
    for strategy in "${strategies[@]}"
    do
        for dts in "${viral[@]}"
        do
            seed=$(make_seed $dts $rep)

            echo "Fine-tuning ESM-2 Viral 650M on $dts with strategy $strategy and seed $rep"
            
            CUDA_VISIBLE_DEVICES=2 python scripts/LitESM2ne_partial_trainer.py \
                -i "data/${source}/metadata/${dts}.csv" \
                -o "experiments/fineTune/${checkpoint}/${source}/${strategy}/" \
                --seed $seed \
                --split_strategy $strategy \
                --checkpoint "$checkpoint"

            echo ""
        done
    done
done

###################################################################
