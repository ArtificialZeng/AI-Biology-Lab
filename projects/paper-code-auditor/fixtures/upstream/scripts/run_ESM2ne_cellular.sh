#!/bin/bash
set -e

# usage
#bash run_ESM2ne_nonviral.sh 

##############################################################################################################
#                                      Cellular sequences
#                                    Riesselman 2018, singles. 
##############################################################################################################
# Apart from BRCA1 batasets (batchsize of 16) I used batchsize of 32 to all others.
cellular=(
    'AMIE_PSEAE_Whitehead'            'DLG4_RAT_Ranganathan2012'            'RL401_YEAST_Bolon2014'
    'B3VI55_LIPSTSTABLE'              'GAL4_YEAST_Shendure2015'             'RL401_YEAST_Fraser2016'
    'B3VI55_LIPST_Whitehead2015'      'HSP82_YEAST_Bolon2016'               'SUMO1_HUMAN_Roth2017'
    'BG_STRSQ_hmmerbit'               'IF1_ECOLI'                           'TIM_SULSO'
    'BLAT_ECOLX_Ostermeier2014'       'KKA2_KLEPN_Mikkelsen2014'            'TIM_THEMA'
    'BLAT_ECOLX_Palzkill2012'         'MK01_HUMAN_Johannessen'              'TIM_THETH'
    'BLAT_ECOLX_Ranganathan2015'      'MTH3_HAEAESTABILIZED_Tawfik2015'     'TPK1_HUMAN_Roth2017'
    'BLAT_ECOLX_Tenaillon2013'        'PABP_YEAST_Fields2013_singles'       'TPMT_HUMAN_Fowler2018'
    'BRCA1_HUMAN_BRCT'                'PTEN_HUMAN_Fowler2018'               'UBC9_HUMAN_Roth2017'
    'BRCA1_HUMAN_RING'                'RASH_HUMAN_Kuriyan'                  'UBE4B_MOUSE_Klevit2013_singles'
    'CALM1_HUMAN_Roth2017'            'RL401_YEAST_Bolon2013'               'YAP1_HUMAN_Fields2012_singles'         
    )
##############################################################################################################

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
source="cellular"

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

##############################################################################################################