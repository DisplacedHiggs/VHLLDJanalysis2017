# LLDJstandalones
ntuple based analysis package for long lived displaced jet analyses

## Download

```bash
# Fermilab uses tcsh by default even though it has bash! 
# This framework is based in bash and 
# technically maybe you don't need this,
# but tcshers be warned
bash --login

# Set up the area
export SCRAM_ARCH=slc6_amd64_gcc630;
scram pro -n LLDJ_slc6_630_CMSSW_9_4_9_cand2 CMSSW CMSSW_9_4_9_cand2;
cd LLDJ_slc6_630_CMSSW_9_4_9_cand2/src;
cmsenv;

## CMSSW imports and customizations
git cms-init --upstream-only
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
scramv1 build -j 20
git cms-merge-topic cms-met:METFixEE2017_949
scramv1 build -j 20

## LLDJstandalones Framework checkout

# first fork the repository to make your own workspace
git clone https://github.com/tmrhombus/LLDJntuples2017.git

 # If you want to check out a specific branch
 # git fetch origin
 # git branch -v -a # list branches available, find yours
 # git checkout -b NAMEOFBRANCH origin/NAMEOFBRANCH 
