import Bio 
from Bio.Data import CodonTable
from Bio.Seq import Seq 
from Bio.Seq import MutableSeq 
from Bio.SeqUtils import GC
import lzma

#SIRT4 = Seq("ATGCCAAATG") // this is responcible for the oxidation of fat mol.(Not needed at this time)
#lzma.compress // for data compression. (may be useful given sequence length)
  
SIRT1_seq = Seq(""" gcatctcctc ctccctctcc ccgggctcct actggcctga ggttgagggc ggctgggggc
       61 tcggggcagg ctccgcggcg ttcccctccc caccccggcc ctccgttcag ccgcgctcct
      121 ccggggctgc ggttcctact gcgcgagctg ccagtggatt cgctcttttc ctccgtccgt
      181 ggcccgcctg ggcggccttg ttctttccgc agcagccaga taaccttctg ttcggtgatg
      241 aaattatcac taatggtttt cattcctgtg aaagtgatga ggaggataga gcctcacatg
      301 caagctctag tgactggact ccaaggccac ggataggtgt ctgtttcatg tggaatacct
      361 gacttcaggt caagggatgg tatttatgct cgccttgctg tagacttccc agatcttcca
      421 gatcctcaag cgatgtttga tattgaatat ttcagaaaag atccaagacc attcttcaag
      481 tttgcaaagg aaatatatcc tggacaattc cagccatctc tctgtcacaa attcatagcc
      541 ttgtcagata aggaaggaaa actacttcgc aactataccc agaacataga cacgctggaa
      601 caggttgcgg gaatccaaag gataattcag tgtcatggtt cctttgcaac agcatcttgc
      661 ctgatttgta aatacaaagt tgactgtgaa gctgtacgag gagatatttt taatcaggta
      721 gttcctcgat gtcctaggtg cccagctgat gaaccgcttg ctatcatgaa accagagatt
      781 gtgttttttg gtgaaaattt accagaacag tttcatagag ccatgaagta tgacaaagat
      841 gaagttgacc tcctcattgt tattgggtct tccctcaaag taagaccagt agcactaatt
      901 ccaagttcca taccccatga agtgcctcag atattaatta atagagaacc tttgcctcat
      961 ctgcattttg atgtagagct tcttggagac tgtgatgtca taattaatga attgtgtcat
     1021 aggttaggtg gtgaatatgc caaactttgc tgtaaccctg taaagctttc agaaattact
     1081 gaaaaacctc cacgaacaca aaaagaattg gcttatttgt cagagttgcc acccacacct
     1141 cttcatgttt cagaagactc aagttcacca gaaagaactt caccaccaga ttcttcagtg
     1201 attgtcacac ttttagacca agcagctaag agtaatgatg atttagatgt gtctgaatca
     1261 aaaggttgta tggaagaaaa accacaggaa gtacaaactt ctaggaatgt tgaaagtatt
     1321 gctgaacaga tggaaaatcc ggatttgaag aatgttggtt ctagtactgg ggagaaaaat
     1381 gaaagaactt cagtggctgg aacagtgaga aaatgctggc ctaatagagt ggcaaaggag
     1441 cagattagta ggcggcttga tggtaatcag tatctgtttt tgccaccaaa tcgttacatt
     1501 ttccatggcg ctgaggtata ttcagactct gaagatgacg tcttatcctc tagttcttgt
     1561 ggcagtaaca gtgatagtgg gacatgccag agtccaagtt tagaagaacc catggaggat
     1621 gaaagtgaaa ttgaagaatt ctacaatggc ttagaagatg agcctgatgt tccagagaga
     1681 gctggaggag ctggatttgg gactgatgga gatgatcaag aggcaattaa tgaagctata
     1741 tctgtgaaac aggaagtaac agacatgaac tatccatcaa acaaatcata gtgtaataat
     1801 tgtgcaggta caggaattgt tccaccagca ttaggaactt tagcatgtca aaatgaatgt
     1861 ttacttgtga actcgataga gcaaggaaac cagaaaggtg taatatttat aggttggtaa
     1921 aatagattgt ttttcatgga taatttttaa cttcattatt tctgtacttg tacaaactca
     1981 acactaactt tttttttttt aaaaaaaaaa aggtactaag tatcttcaat cagctgttgg
     2041 tcaagactaa ctttctttta aaggttcatt tgtatgataa attcatatgt gtatatataa
     2101 ttttttttgt tttgtctagt gagtttcaac atttttaaag ttttcaaaaa gccatcggaa
     2161 tgttaaatta atgtaaaggg aacagctaat ctagaccaaa gaatggtatt ttcacttttc
     2221 tttgtaacat tgaatggttt gaagtactca aaatctgtta cgctaaactt ttgattcttt
     2281 aacacaatta tttttaaaca ctggcatttt ccaaaactgt ggcagctaac tttttaaaat
     2341 ctcaaatgac atgcagtgtg agtagaagga agtcaacaat atgtggggag agcactcggt
     2401 tgtctttact tttaaaagta atacttggtg ctaagaattt caggattatt gtatttacgt
     2461 tcaaatgaag atggcttttg tacttcctgt ggacatgtag caatgtctat attggctcat
     2521 aaaactaacc tgaaaaacaa ataaatgctt tggaaatgtt tcagttgctt tagaaacatt
     2581 agtgcctgcc tggatcccct tagttttgaa atatttgcca ttgttgttta aatacctatc
     2641 actgtggtag agcttgcatt gatcttttcc acaagtatta aactgccaaa atgtgaatat
     2701 gcaaagcctt tctgaatcta taataatggt acttctactg gggagagtgt aatattttgg
     2761 actgctgttt tccattaatg aggagagcaa caggcccctg attatacagt tccaaagtaa
     2821 taagatgtta attgtaattc agccagaaag tacatgtctc ccattgggag gatttggtgt
     2881 taaataccaa actgctagcc ctagtattat ggagatgaac atgatgatgt aacttgtaat
     2941 agcagaatag ttaatgaatg aaactagttc ttataattta tctttattta aaagcttagc
     3001 ctgccttaaa actagagatc aactttctca gctgcaaaag cttctagtct ttcaagaagt
     3061 tcatacttta tgaaattgca cagtaagcat ttatttttca gaccattttt gaacatcact
     3121 cctaaattaa taaagtattc ctctgttgct ttagtattta ttacaataaa aagggtttga
     3181 aatatagctg ttctttatgc ataaaacacc cagctaggac cattactgcc agagaaaaaa
     3241 atcgtattga atggccattt ccctacttat aagatgtctc aatctgaatt tatttggcta
     3301 cactaaagaa tgcagtatat ttagttttcc atttgcatga tgtttgtgtg ctatagatga
     3361 tattttaaat tgaaaagttt gttttaaatt atttttacag tgaagactgt tttcagctct
     3421 ttttatattg tacatagtct tttatgtaat ttactggcat atgttttgta gactgtttaa
     3481 tgactggata tcttccttca acttttgaaa tacaaaacca gtgtttttta cttgtacact
     3541 gttttaaagt ctattaaaat tgtcatttga cttttttctg ttaactta """)  
SIRT6_seq = Seq(""" attgttcccg tggggcagtc gaggatgtcg gtgaattacg cggcggggct gtcgccgtac
       61 gcggacaagg gcaagtgcgg cctcccggag atcttcgacc ccccggagga gctggagcgg
      121 aaggtgtggg aactggcgag gctggtctgg cagtcttcca gtgtggtgtt ccacacgggt
      181 gccggcatca gcactgcctc tggcatcccc gacttcaggg gtccccacgg agtctggacc
      241 atggaggagc gaggtctggc ccccaagttc gacaccacct ttgagagcgc gcggcccacg
      301 cagacccaca tggcgctggt gcagctggag cgcgtgggcc tcctccgctt cctggtcagc
      361 cagaacgtgg acgggctcca tgtgcgctca ggcttcccca gggacaaact ggcagagctc
      421 cacgggaaca tgtttgtgga agaatgtgcc aagtgtaaga cgcagtacgt ccgagacaca
      481 gtcgtgggca ccatgggcct gaaggccacg ggccggctct gcaccgtggc taaggcaagg
      541 gggctgcgag cctgcaggaa cgccgacctg tccatcacgc tgggtacatc gctgcagatc
      601 cggcccagcg ggaacctgcc gctggctacc aagcgccggg gaggccgcct ggtcatcgtc
      661 aacctgcagc ccaccaagca cgaccgccat gctgacctcc gcatccatgg ctacgttgac
      721 gaggtcatga cccggctcat gaagcacctg gggctggaga tccccgcctg ggacggcccc
      781 cgtgtgctgg agagggcgct gccacccctg ccccgcccgc ccacccccaa gctggagccc
      841 aaggaggaat ctcccacccg gatcaacggc tctatccccg ccggccccaa gcaggagccc
      901 tgcgcccagc acaacggctc agagcccgcc agccccaaac gggagcggcc caccagccct
      961 gccccccaca gaccccccaa aagggtgaag gccaaggcgg tccccagctg accagggtgc
     1021 ttggggaggg tggggctttt tgtagaaact gtggattctt tttctctcgt ggtctcactt
     1081 tgttacttgt ttctgtcccc gggagcctca gggctctgag agctgtgctc caggccaggg
     1141 gttacacctg ccctccgtgg tccctccctg ggctccaggg gcctctggtg cggttccggg
     1201 aagaagccac accccagagg tgacaggtga gcccctgcca caccccagcc tctgacttgc
     1261 tgtgttgtcc agaggtgagg ctgggccctc cctggtctcc agcttaaaca ggagtgaact
     1321 ccctctgtcc ccagggcctc ccttctgggc cccctacagc ccaccctacc cctcctccat
     1381 gggccctgca ggaggggaga cccaccttga agtgggggat cagtagaggc ttgcactgcc
     1441 tttggggctg gagggagacg tgggtccacc aggcttctgg aaaagtcctc aatgcaataa
     1501 aaacaatttc tttcttgca """)
SIRT7_seq = Seq(""" 1 ctgccgtgtg aggcggaagc ggaagagcag gtctccaggg gagcgatggc agccgggggt
       61 ctgagccgct ccgagcgcaa agcggcggag cgggtccgga ggttgcggga ggagcagcag
      121 agggagcgcc tccgccaggt gtcgcgcatc ctgaggaagg cggcggcgga gcgcagcgcc
      181 gaggagggcc ggctgctggc cgagagcgcg gacctggtaa cggagctgca gggccggagc
      241 cggcggcgcg agggcctgaa gcggcggcag gaggaggtgt gcgacgaccc ggaggagctg
      301 cgggggaagg tccgggagct ggccagcgcc gtccggaacg ccaaatactt ggtcgtctac
      361 acaggcgcgg gaatcagcac ggcagcgtct atcccagact accggggccc taatggagtg
      421 tggacactgc ttcagaaagg gagaagcgtt agtgctgccg acctgagcga ggccgagcca
      481 accctcaccc acatgagcat cacccgtctg catgagcaga agctggtgca gcatgtggtg
      541 tctcagaact gtgacgggct ccacctgagg agtgggctgc cgcgcacggc catctccgag
      601 ctccacggga acatgtacat tgaagtctgt acctcctgcg ttcccaacag ggagtacgtg
      661 cgggtgttcg atgtgacgga gcgcactgcc ctccacagac accagacagg ccggacctgc
      721 cacaagtgtg ggacccagct gcgggacacc attgtgcact ttggggagag ggggacgttg
      781 gggcagcctt tgaactggga agcggcgacc gaggctgcca gcagagcaga caccatcctg
      841 tgtctagggt ccagcctgaa ggttctaaag aagtacccac gcctctggtg catgaccaag
      901 ccccctagcc ggcggccgaa gctttacatc gtgaacctgc agtggacccc gaaggatgac
      961 tgggctgccc tgaagctaca tgggaagtgt gatgacgtca tgcggctcct catggccgag
     1021 ctgggcttgg agatccccgc ctatagcagg tggcaggatc ccattttctc actggcgact
     1081 cccctgcgtg ctggtgaaga aggcagccac agtcggaagt cgctgtgcag aagcagagag
     1141 gaggccccgc ctggggaccg gggtgcaccg cttagctcgg cccccatcct agggggctgg
     1201 tttggcaggg gctgcacaaa acgcacaaaa aggaagaaag tgacgtaatc acgtgctcga
     1261 tgaagaacag ttggcacttt gcagatggcc agtgtcacgg tgaaggctgg gttgccccca
     1321 cgggtctagg gagaacgaac tctttgggga tgacattttc accgtgacat ttttagccat
     1381 ttgtccttga ggaagcccct tgcactgctg cggttgtacc ctgatacggc ctggccatcg
     1441 aggacacctg cccatccggc ctctgtgtca agaggtggca gccgcacctt tctgtgagaa
     1501 cggaactcgg gttatttcag ccccggcctg cagagtggaa gcgcccagcg gcctttcctc
     1561 gctcaccagg ccagtctcag ggcctcaccg tatttctact actacttaat gaaaaagtgt
     1621 gaactttata gaatcctctc tgtactggat gtgcggcaga ggggtggctc cgagcctcgg
     1681 ctctatgcag acctttttat ttctattaaa cgtttctgca ctggc """)

#SIRT1 promotes homologus recombination in cells + recombo in DNA breaks  
#SIRT6 Chromatin-associated protein, base repair of DNA
#SIRT7 employed to repair double strand breaks
 
transc_rna1 = SIRT1_seq.transcribe() 
transc_rna6 = SIRT6_seq.transcribe()
transc_rna7 = SIRT7_seq.transcribe()
#print(transc_rna1, transc_rna6, transc_rna7) 

