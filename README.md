# Benchmark of software tools for prokaryotic chromosomal interaction domain identification

Mikhail D Magnitov, Veronika S Kuznetsova, Sergey V Ulianov, Sergey V Razin, Alexander V Tyakht (2020). Benchmark of software tools for prokaryotic chromosomal interaction domain identification.

doi: [10.1093/bioinformatics/btaa555](https://doi.org/10.1093/bioinformatics/btaa555)

**Datasets used:**
* *Caulobacter crescentus* (Le *et al.*, *Science*, 2013)
* *Bacillus subtilis* (Wang *et al.*, *Genes Dev.*, 2015; Marbouty *et al.*, *Mol. Cell*, 2015)
* *Escherichia coli* (Lioy *et al.*, *Cell*, 2018)
* *Mycoplasma pneumoniae* (Trussart *et al.*, *Nat. Comm.*, 2017)
* *Sulfolobus acidocaldarius* (Takemata *et al.*, *Cell*, 2019)

**Compared domain-calling tools:**
* [Armatus](https://github.com/kingsfordgroup/armatus/) (Filippova *et al.*, *Algorithms Mol Biol.*, 2014)
* [Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead) (Durand *et al.*, *Cell Systems*, 2016)
* [CaTCH](https://github.com/zhanyinx/CaTCH_R) (Zhan *et al.*, *Genome Res.*, 2017)
* [CHDF](https://link.springer.com/article/10.1007/s40484-015-0047-9) (Wang *et al.*, *Quant Biol*, 2015)
* [chromoR](https://cran.r-project.org/web/packages/chromoR/index.html) (Shavit and Lio, *Mol. Biosyst.*, 2014)
* [Chromosight](https://github.com/koszullab/chromosight) (Matthey-Doret *et al.*, *bioRxiv*, 2020)
* [ClusterTAD](https://github.com/BDM-Lab/ClusterTAD) (Oluwadare and Cheng, *BMC Bioinformatics*, 2017)
* [deDoc](https://github.com/yinxc/structural-information-minimisation) (Li *et al.*, *Nat. Comm.*, 2018)
* Directionality Index (Dixon *et al.*, *Nature*, 2012; Le *et al.*, *Science*, 2013)
* [EAST](https://github.com/ucrbioinfo/EAST) (Roayaei Ardakany and Lonardi, *17th International Workshop on Algorithms in Bioinformatics*, 2017)
* [GMAP](https://github.com/wbaopaul/rGMAP) (Yu *et al.*, *Nat. Comm.*, 2017)
* [HiCExplorer](https://github.com/deeptools/HiCExplorer) (Ramirez *et al.*, *Nat. Comm.*, 2018)
* [HiCseg](https://CRAN.R-project.org/package=HiCseg) (Levy-Leduc *et al.*, *Bioinformatics*, 2014)
* [IC-Finder](https://github.com/bcm-uga/IC-Finder) (Haddad *et al.*, *Nucleic Acids Res.*, 2017)
* [Insulation Score](https://github.com/dekkerlab/crane-nature-2015/) (Crane *et al.*, *Nature*, 2015)
* [Lavaburst](https://github.com/nvictus/lavaburst)
* [MrTADFinder](https://github.com/gersteinlab/MrTADFinder) (Yan *et al.*, *PLoS Comput. Biol.*, 2017)
* [OnTAD](https://github.com/anlin00007/OnTAD) (An *et al.*, *Genome Biol.*, 2019)
* [spectral](https://github.com/laseaman/4D_Nucleome_Analysis_Toolbox) (Chen *et al.*, *Bioinformatics*, 2016)
* [SpectralTAD](https://github.com/dozmorovlab/SpectralTAD) (Cresswell *et al.*, *bioRxiv*, 2019)
* [TADbit](https://github.com/3DGenomes/TADbit) (Serra *et al.*, *PLoS Comput. Biol.*, 2017)
* [TADpole](https://github.com/3DGenomes/TADpole) (Soler-Villa *et al.*, *Nucleic Acids Res.*, 2020)
* [TADtool](https://github.com/vaquerizaslab/tadtool) (Kruse *et al.*, *Bioinformatics*, 2016)
* [TADtree](http://compbio.cs.brown.edu/projects/tadtree/) (Weinreb *et al.*, *Bioinformatics*, 2015)
* [TopDom](https://github.com/HenrikBengtsson/TopDom) (Shin *et al.*, *Nucleic Acids Res.*, 2016)

**Available code and data:**
* Hi-C and RNA-seq data processing pipelines
* Scripts to run domain callers used for benchmarking
* Annotated domain positions used in the analysis
* Scripts to analyse and validate the annotated domains
* Scripts to digest genomes and simulate Hi-C read pairs
