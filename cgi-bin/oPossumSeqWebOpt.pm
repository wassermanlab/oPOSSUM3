#
# Sequence-based oPOSSUM specific options
#

use constant DFLT_CORE_MIN_IC   => 8;

#use constant DFLT_BACKGROUND    => 'MM0444_min05_peak_seqs.fa';
use constant DFLT_BACKGROUND    => 'MM0444_min05_5000_rand_peak_seqs.fa';

use constant BG_SEQ_SET_KEYS   => [
    'mmFibro', 'mmLiver', 'mmES', 'mmMarrow', 'mmMixed'
];

use constant BG_SEQ_SET_NAMES  => {
    mmFibro     => 'Mouse fibroblast cells',
    mmLiver     => 'Mouse liver cells',
    mmES        => 'Mouse ES cells',
    mmMarrow    => 'Mouse bone marrow cells',
    mmMixed     => 'Mouse combined cell lines'
};

use constant BG_SEQ_SET_FILES  => {
    mmFibro     => 'fibroblast_rand5000_seq.fa',
    mmLiver     => 'liver_rand5000_seq.fa',
    mmES        => 'ES_rand5000_seq.fa',
    mmMarrow    => 'bonemarrow_rand5000_seq.fa',
    mmMixed     => 'mixture_ES-Liver-Bonemarrow-eFibroblast_rand5000_seq.fa'
};

# TF DB Information
use constant JASPAR_COLLECTIONS => [
    'PHYLOFACTS',
    'FAM',
    'PBM',
    'PBM_HOMEO',
    'PBM_HLH'
];

use constant JASPAR_CORE_TAX_GROUPS => [
    'vertebrates',
    #'urochordates',
    'insects',
    'nematodes',
    'plants',
    'fungi'
];

# oPOSSUM default parameters (for custom analysis).
use constant MIN_THRESHOLD      => 70;

1;
