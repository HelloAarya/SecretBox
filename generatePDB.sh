#!/bin/bash

function show_help () {
cat <<-END

    Pipeline for local trRosetta model generation and validation
    
    Author: Aarya Venkat and Wayland Yeung

    The following program takes in a fasta file and converts it into an a3m file.
    It then uses this a3m alignment to produce an npz array of dist/angles. Finally, 
    it uses the npz array and the fasta sequence you'd like to predict a model from
    generate multiple models of your sequence. As a note: be careful continuing if 
    the contact score if below 0.4 or if the fasta seq does not match the NPZ file 
    length (check.py)

    Requirements:
    # Python3
    #     sudo apt install python3
    # TensorFlow 1.14
    #     sudo python3 -m pip install tensorflow=1.14
    # Numpy
    #     sudo python3 -m pip install numpy
    # PyRosetta
    #     cd pyrosetta/setup && sudo python3 setup.py install

    Usage:
    ./generatePDB.sh <options> fasta.fa

    -h show help and exit
        toggle
    -o output name
        argument
    -r remove inserts
        toggle
    -d database (hmm)
        argument
    -e expand the alignment
        toggle
    -a alignment file (fasta)
        argument
    -n number of models
        argument
    -t number of threads
        argument

END
}

if [[ $# -eq 0 ]]; then

    ./generatePDB.sh -h;

    exit 0
fi


output_name=''
has_output_name=false

remove_inserts=false

database=''
has_database=false

expand_alignment=false

alignment_file=''
has_alignment_file=false

num_models=10
num_threads=1

while getopts ho:rd:ea:n:t: opt; do
    case $opt in
        h) show_help; exit 1                               ;;
        o) output_name=$OPTARG; has_output_name=true       ;;
        r) remove_inserts=true                             ;;
        d) database=$OPTARG; has_database=true             ;;
        e) expand_alignment=true                           ;;
        a) alignment_file=$OPTARG; has_alignment_file=true ;;
        n) num_models=$OPTARG                              ;;
        t) num_threads=$OPTARG                             ;;
    esac
done

shift $(( OPTIND - 1 ))
input_file=$1

################################################################################
################################## create the msa

if $has_database; then
  true # do nothing
else
  database=/auto/share/db/pdb70/pdb70
fi

if $has_output_name; then
  true # do nothing
else
  output_name=trr_${input_file}
  # mkdir -r `dirname $output_name`
fi

prod_input=${output_name}.input
prod_aln=${output_name}.align

./scripts/reformat.pl fas a3m $input_file $prod_input # standardized input file

if $has_alignment_file; then
  ./scripts/reformat.pl fas a3m $alignment_file $prod_aln
  if $expand_alignment; then
    ./hhsuite/bin/hhblits -i $prod_aln -oa3m ${prod_aln}.expand -d $database
    prod_aln=${prod_aln}.expand
  fi
else
  ./hhsuite/bin/hhblits -i $prod_input -oa3m $prod_aln -d $database
fi

if $remove_inserts; then
  ./scripts/reformat.pl a3m a3m $prod_aln ${prod_aln}.noins -r
  prod_aln=${prod_aln}.noins
fi

./scripts/rm_newline.pl $prod_aln > ${prod_aln}.temp
mv ${prod_aln}.temp $prod_aln

################################################################################
################################## create the npz

prod_model=${output_name}.npz
prod_report=${output_name}.log

python3 network/predict.py -m model2019_07 $prod_aln $prod_model

# Calculates probability of Top Contacts from npz file (TMscore of the best structure).
python3 scripts/top_prob.py $prod_model > $prod_report


# Evaluates the sequence length of the fasta vs npz file
python3 scripts/check.py $prod_model $prod_input >> $prod_report

################################################################################
################################## generate models

echo -e '\n\n\n########################## MODEL SCORES\n'>> $prod_report

for i in $(seq 1 $num_models); do
  current_model=${output_name}_$i.pdb
  echo "python3 scripts/trRosetta.py $prod_model $prod_input $current_model; printf "$current_model " >> $prod_report; grep pose $current_model | cut -d ' ' -f 24 >> $prod_report"
done | parallel -uj $num_threads
