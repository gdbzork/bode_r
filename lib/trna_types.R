# Some basic types for tRNA manipulations

species <- c('cfa', 'hep', 'hsa', 'mdo', 'mml', 'mmu', 'rno')

isotypes <- c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
              'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'SeC', 'Ser', 'Thr',
              'Trp', 'Tyr', 'Val')

anticodons.all <- c('Ala_AGC', 'Ala_CGC', 'Ala_TGC', 'Arg_ACG', 'Arg_CCG',
                    'Arg_CCT', 'Arg_TCG', 'Arg_TCT', 'Asn_GTT', 'Asp_GTC',
                    'Cys_GCA', 'Gln_CTG', 'Gln_TTG', 'Glu_CTC', 'Glu_TTC',
                    'Gly_ACC', 'Gly_CCC', 'Gly_GCC', 'Gly_TCC', 'His_GTG',
                    'Ile_AAT', 'Ile_GAT', 'Ile_TAT', 'Leu_AAG', 'Leu_CAA',
                    'Leu_CAG', 'Leu_GAG', 'Leu_TAA', 'Leu_TAG', 'Lys_CTT',
                    'Lys_TTT', 'Met_CAT', 'Phe_GAA', 'Pro_AGG', 'Pro_CGG',
                    'Pro_GGG', 'Pro_TGG', 'SeC_TCA', 'Ser_AGA', 'Ser_CGA',
                    'Ser_GCT', 'Ser_GGA', 'Ser_TGA', 'Thr_AGT', 'Thr_CGT',
                    'Thr_TGT', 'Trp_CCA', 'Tyr_GTA', 'Val_AAC', 'Val_CAC',
                    'Val_TAC')

anticodons.core <- c('Ala_AGC', 'Ala_CGC', 'Ala_TGC', 'Arg_ACG', 'Arg_CCG',
                     'Arg_CCT', 'Arg_TCG', 'Arg_TCT', 'Asn_GTT', 'Asp_GTC',
                     'Cys_GCA', 'Gln_CTG', 'Gln_TTG', 'Glu_CTC', 'Glu_TTC',
                     'Gly_CCC', 'Gly_GCC', 'Gly_TCC', 'His_GTG', 'Ile_AAT',
                     'Ile_TAT', 'Leu_AAG', 'Leu_CAA', 'Leu_CAG', 'Leu_TAA',
                     'Leu_TAG', 'Lys_CTT', 'Lys_TTT', 'Met_CAT', 'Phe_GAA',
                     'Pro_AGG', 'Pro_CGG', 'Pro_TGG', 'SeC_TCA', 'Ser_AGA',
                     'Ser_CGA', 'Ser_GCT', 'Ser_TGA', 'Thr_AGT', 'Thr_CGT',
                     'Thr_TGT', 'Trp_CCA', 'Tyr_GTA', 'Val_AAC', 'Val_CAC',
                     'Val_TAC')
