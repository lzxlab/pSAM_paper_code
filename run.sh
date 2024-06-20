##squeue -o "%i %.30j %.8u %.2t %.8M %.6D %.5C %.5m %.15R %Z"

##################model construction##################
#######build by nucleus
##train-valid-test
/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./0.10.train.valid.test.R Nucleus
##test different models and parameters
awk '{print "sbatch 1.0.gridsearch.sh "$1""}' models.txt | bash


##rocfile for each model, train-valid-test, cvs
awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./1.1.train.valid.82.R Nucleus "$1""}' nums.txt | bash
awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./1.1.train.valid.82.R Cytosol "$1""}' nums.txt | bash
awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./1.1.train.valid.82.R Membrane "$1""}' nums.txt | bash
awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./1.1.train.valid.82.R Mitochondrion "$1""}' nums.txt | bash
awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./1.1.train.valid.82.R Secreted "$1""}' nums.txt | bash
awk '{print "sbatch 1.1.cv.sh "$1" "$2""}' locs.nums.txt | bash

awk '{print "/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./0.10.train.valid.test.R "$1""}' locs.txt | bash
awk '{print "sbatch 2.0.tr.va.te.sh "$1""}' locs.txt | bash

##define final model for all locations
awk '{print "cat ../0.data/1.dataset.subloc."$1".txt.train ../0.data/1.dataset.subloc."$1".txt.valid ../0.data/1.dataset.subloc."$1".txt.test > ../0.data/1.dataset.subloc."$1""}' locs.txt | bash
awk '{print "sbatch 2.1.alldata.sh "$1""}' locs.txt | bash

##################point mutation and paired mutation##################
#######point mutation for each location
awk '{print "sbatch 4.0.hs.score.sh "$1""}' locs.txt | bash
rm ../4.panMut/merged.mut.maf.score*
awk '{print "sbatch 4.1.predictmut.sh "$1" "$2""}' locs.mut.txt | bash

#######paired mutation for each location
rm ../5.panpairMut/merged.mut.pair.maf.score*
awk '{print "sbatch 4.2.predictpairmut.sh  "$1" "$2""}' locs.pairmut.txt | bash

#######point mutation for each location (just for Nucleus, for mistake)
sh 4.0.hs.score.sh Nucleus
rm ../4.panMut/merged.mut.maf.score.Nucleus
awk '{print "sh 4.1.predictmut.sh "$1" "$2""}' locs.mut.txt | grep Nucleus | bash
#######paired mutation for each location
rm ../5.panpairMut/merged.mut.pair.maf.score.Nucleus
awk '{print "sh 4.2.predictpairmut.sh  "$1" "$2""}' locs.pairmut.txt | grep Nucleus |bash


##################other information##################
#######predict 1000 sequence with existing tools#######
awk '{print "sh 3.0.predict1000.sh "$1""}' locs.txt | grep -v Nucleus | bash

#######predict LPT for all organisms#######
##cutoff set
##Nucleus:0.3149788, Cytosol:0.6039815, Membrane:0.0868726, 
##Mitochondrion:0.1380336, Secreted:0.05454149, 
#awk '{print "sh 7.0.predict.allpros.sh "$1" Nucleus 0.3149788"}' ../7.eachpro/inuloc.dataset.scores | xargs -iCMD -P10 bash -c CMD
sbatch 7.0.predict.allpros.sbatch "$1" Cytosol 0.6039815
sbatch 7.0.predict.allpros.sbatch "$1" Membrane 0.0868726
sbatch 7.0.predict.allpros.sbatch "$1" Mitochondrion 0.1380336
sbatch 7.0.predict.allpros.sbatch "$1" Secreted 0.05454149

cat ../7.eachpro/Nucleus/*.delta.scores > ../7.eachpro/hs.CPT.delta.scores.Nucleus
cat ../7.eachpro/Cytosol/*.delta.scores > ../7.eachpro/hs.CPT.delta.scores.Cytosol
cat ../7.eachpro/Membrane/*.delta.scores > ../7.eachpro/hs.CPT.delta.scores.Membrane
cat ../7.eachpro/Mitochondrion/*.delta.scores > ../7.eachpro/hs.CPT.delta.scores.Mitochondrion
cat ../7.eachpro/Secreted/*.delta.scores > ../7.eachpro/hs.CPT.delta.scores.Secreted

cat ../7.eachpro/Nucleus/*.delta.scores.regions > ../7.eachpro/hs.CPT.delta.scores.regions.Nucleus
cat ../7.eachpro/Cytosol/*.delta.scores.regions > ../7.eachpro/hs.CPT.delta.scores.regions.Cytosol
cat ../7.eachpro/Membrane/*.delta.scores.regions > ../7.eachpro/hs.CPT.delta.scores.regions.Membrane
cat ../7.eachpro/Mitochondrion/*.delta.scores.regions > ../7.eachpro/hs.CPT.delta.scores.regions.Mitochondrion
cat ../7.eachpro/Secreted/*.delta.scores.regions > ../7.eachpro/hs.CPT.delta.scores.regions.Secreted

##################get regions for each locs##################
python 7.1.get.regions.py hs.CPT.delta.scores.regions.sigdelta.Cytosol
python 7.1.get.regions.py hs.CPT.delta.scores.regions.sigdelta.Membrane
python 7.1.get.regions.py hs.CPT.delta.scores.regions.sigdelta.Mitochondrion
python 7.1.get.regions.py hs.CPT.delta.scores.regions.sigdelta.Secreted


##################NMC terminal trunc##################
python 6.0.predict.trunc.py Cytosol
python 6.0.predict.trunc.py Membrane
python 6.0.predict.trunc.py Mitochondrion
python 6.0.predict.trunc.py Secreted

##################regions truncated##################
python 6.1.predict.regions.py mito.human.pos.txt Mitochondrion
python 6.1.predict.regions.py sp.human.pos.txt Secreted

##################sig regions truncated##################
python 6.1.predict.regions.py hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions Cytosol
python 6.1.predict.regions.py hs.CPT.delta.scores.regions.sigdelta.Membrane.regions Membrane
python 6.1.predict.regions.py hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions Mitochondrion
python 6.1.predict.regions.py hs.CPT.delta.scores.regions.sigdelta.Secreted.regions Secreted

##################regions seq retrive##################
python 7.2.get.fasta.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions
python 7.2.get.fasta.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions
python 7.2.get.fasta.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions
python 7.2.get.fasta.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions

##################motif enrichment##################
seqkit seq -g -m 5 ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta > ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5
seqkit seq -g -m 5 ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.fasta > ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.fasta5
seqkit seq -g -m 5 ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.fasta > ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.fasta5
seqkit seq -g -m 5 ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.fasta > ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.fasta5
seqkit seq -g -m 5 ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.fasta > ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.fasta5

meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.fasta5 -oc ../7.eachpro/meme510Cytosol -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 10 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.fasta5 -oc ../7.eachpro/meme515Cytosol -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 15 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.fasta5 -oc ../7.eachpro/meme510Membrane -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 10 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.fasta5 -oc ../7.eachpro/meme515Membrane -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 15 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.fasta5 -oc ../7.eachpro/meme510Mitochondrion -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 10 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.fasta5 -oc ../7.eachpro/meme515Mitochondrion -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 15 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.fasta5 -oc ../7.eachpro/meme510Secreted -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 10 -p 3 -V
meme ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.fasta5 -oc ../7.eachpro/meme515Secreted -seed 123456 -objfun de -neg ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta5 -protein -nmotifs 10 -minw 5 -maxw 15 -p 3 -V

##################NLS pSTY Functional score##################
python 7.4.find.position.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.nls.regions ../7.eachpro/pSTY.NLS.position.csv
python 7.4.find.position.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions ../7.eachpro/pSTY.Cytosol.position.csv
python 7.4.find.position.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions ../7.eachpro/pSTY.Membrane.position.csv
python 7.4.find.position.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions ../7.eachpro/pSTY.Mitochondrion.position.csv
python 7.4.find.position.py ../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions ../7.eachpro/pSTY.Secreted.position.csv

##################selected gene mut NLPT##################
rm /home/yukai/projects/CBAmodel/ProNuclear/7.singlemut/tmp/* ../4.panMut/mutplot/*
python /home/yukai/projects/CBAmodel/ProNuclear/7.singlemut/6.predict_singlemut.py 	Q96EP1 P_255_L
cp -r /home/yukai/projects/CBAmodel/ProNuclear/7.singlemut/tmp/* ../4.panMut/mutplot/

##################NLS pSTY mutation##################
python 4.3.nls.psty.mut.py

##################basic analysis and results##################
Rscript analysis.R





