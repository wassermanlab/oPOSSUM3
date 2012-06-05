#
# Sequence-based oPOSSUM web specific options
#

#use constant DFLT_BACKGROUND    => 'MM0444_min05_peak_seqs.fa';
use constant DFLT_BACKGROUND    => 'MM0444_min05_5000_rand_peak_seqs.fa';

use constant BG_SEQ_SET_KEYS   => [
    'mmFibro2500',
    'mmFibro5000',
    'mmLiver2500',
    'mmLiver5000',
    'mmMarrow2500',
    'mmMarrow5000',
    'mmMixed2500',
    'mmMixed5000'
];

use constant BG_SEQ_SET_NAMES  => {
    mmFibro2500     => 'Mouse fibroblast 2500 seqs (GC=44%)',
    mmFibro5000     => 'Mouse fibroblast 5225 seqs (GC=44%)',
    mmLiver2500     => 'Mouse liver 2500 seqs (GC=51%)',
    mmLiver5000     => 'Mouse liver 5000 seqs (GC=51%)',
    mmMarrow2500    => 'Mouse bone marrow 2500 seqs (GC=46%)',
    mmMarrow5000    => 'Mouse bone marrow 5000 seqs (GC=46%)',
    mmMixed2500     => 'Mouse mixed cell lines 2500 seqs (GC=45%)',
    mmMixed5000     => 'Mouse mixed cell lines 5000 seqs (GC=45%)'
};

use constant BG_SEQ_SET_FILES  => {
    mmFibro2500     => 'mouse_fibroblast_2500seq_44percentGC.fa',
    mmFibro5000     => 'mouse_fibroblast_5225seq_44percentGC.fa',
    mmLiver2500     => 'mouse_liver_2500seq_51percentGC.fa',
    mmLiver5000     => 'mouse_liver_5000seq_51percentGC.fa',
    mmMarrow2500    => 'mouse_bonemarrow_2500seq_46percentGC.fa',
    mmMarrow5000    => 'mouse_bonemarrow_5000seq_46percentGC.fa',
    mmMixed2500     => 'mouse_mixture_EScell-Liver-Bonemarrow-eFibroblast_2500seq_45percentGC.fa',
    mmMixed5000     => 'mouse_mixture_EScell-Liver-Bonemarrow-eFibroblast_5000seq_45percentGC.fa'
};

1;
