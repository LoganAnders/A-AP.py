import Bio 
from Bio.Data import CodonTable
from Bio.Seq import Seq 
from Bio.Seq import MutableSeq 
from Bio.SeqUtils import GC
import lzma

#SIRT4 = Seq("ATGCCAAATG") // this is responcible for the oxidation of fat mol.(Not needed at this time)
#lzma.compress // for data compression. (may be useful given sequence length)
  
SIRT1_seq = Seq("""  mfdieyfrkd prpffkfake iypgqfqpsl chkfialsdk egkllrnytq nidtleqvag
        iqriiqchgs fatasclick ykvdceavrg difnqvvprc prcpadepla imkpeivffg
       enlpeqfhra mkydkdevdl livigsslkv rpvalipssi phevpqilin replphlhfd
       vellgdcdvi inelchrlgg eyaklccnpv klseitekpp rtqkelayls elpptplhvs
       edsssperts ppdssvivtl ldqaaksndd ldvseskgcm eekpqevqts rnvesiaeqm
       enpdlknvgs stgeknerts vagtvrkcwp nrvakeqisr rldgnqylfl ppnryifhga
       evysdseddv lsssscgsns dsgtcqspsl eepmedesei eefynglede pdvperagga
       gfgtdgddqe aineaisvkq evtdmnypsn ks """)  
SIRT6_seq = Seq(""" msvnyaagls pyadkgkcgl peifdppeel erkvwelarl vwqsssvvfh tgagistasg
        ipdfrgphgv wtmeerglap kfdttfesar ptqthmalvq lervgllrfl vsqnvdglhv
        rsgfprdkla elhgnmfvee cakcktqyvr dtvvgtmglk atgrlctvak arglracrna
       dlsitlgtsl qirpsgnlpl atkrrggrlv ivnlqptkhd rhadlrihgy vdevmtrlmk
       hlgleipawd gprvleralp plprpptpkl epkeesptri ngsipagpkq epcaqhngse
       paspkrerpt spaphrppkr vkakavps  """)
SIRT7_seq = Seq(""" maagglsrse rkaaervrrl reeqqrerlr qvsrilrkaa aersaeegrl laesadlvte
        lqgrsrrreg lkrrqeevcd dpeelrgkvr elasavrnak ylvvytgagi staasipdyr
       gpngvwtllq kgrsvsaadl seaeptlthm sitrlheqkl vqhvvsqncd glhlrsglpr
       taiselhgnm yievctscvp nreyvrvfdv tertalhrhq tgrtchkcgt qlrdtivhfg
       ergtlgqpln weaateaasr adtilclgss lkvlkkyprl wcmtkppsrr pklyivnlqw
       tpkddwaalk lhgkcddvmr llmaelglei paysrwqdpi fslatplrag eegshsrksl
       crsreeappg drgaplssap ilggwfgrgc tkrtkrkkvt  """)

analysed_SIRT1 = ProteinAnalysis(SIRT1_seq)
analysed_SIRT6 = ProteinAnalysis(SIRT6_seq)
analysed_SIRT7 = ProteinAnalysis(SIRT7_seq) 

#SIRT1 promotes homologus recombination in cells + recombo in DNA breaks  
#SIRT6 Chromatin-associated protein, base repair of DNA
#SIRT7 employed to repair double strand breaks
 
analysed_SIRT1.molecular_weight() 
analysed_SIRT1.get_amino_acids_percent() 



