#!/bin/bash
trap "exit" INT

WHAT=$1;
if [ "$#" -ne 1 ]; then
    echo "steerTOPJS.sh <SEL/PLOTSEL/WWWSEL>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    exit 1;
fi

export LSB_JOB_REPORT_MAIL=N


queue=longlunch
githash=b312177
lumi=35922
lumiSpecs="" #--lumiSpecs EE:11391"
lumiUnc=0.027
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
run_summaryeosdir=/store/group/cmst3/user/mseidel/analysis/TopJetShapes/b312177
summaryeosdir=/eos/cms/store/group/cmst3/user/mseidel/analysis/TopJetShapes/b312177
wwwdir=/eos/user/${myletter}/${whoami}/www/cms/TopJS/


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        scram b -j 8 && python scripts/runLocalAnalysis.py -i ${eosdir}/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root --tag MC13TeV_TTJets -o analysis.root --era era2016 -m TOPJetShape::RunTopJetShape --debug;
        ;;

    NORMCACHE )
        python scripts/produceNormalizationCache.py -i ${eosdir} -o data/era2016/genweights.root;
        python scripts/produceNormalizationCache.py -i ${eosdir}_qcd -o data/era2016/genweights.root --update;
        ;;

    FULLSEL )
        #ttbar, modeling systematics, background samples
        python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --skip Data13TeV_Double,Data13TeV_MuonEG,TTJets2l2nu,m166,m169,m175,m178,widthx,TTTT --farmappendix samples;
        #QCD samples
        python scripts/runLocalAnalysis.py -i ${eosdir}_qcd -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --farmappendix qcd;
        #Experimental uncertainties
        python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --only MC13TeV_TTJets --systVar all --exactonly --farmappendix expsyst;
        #ttbar GEN samples
        python scripts/runLocalAnalysis.py -i /store/group/cmst3/user/mseidel/ReReco2016/b312177_GEN_merged -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --farmappendix samplesGEN;
        ;;

    FULLSELCENTRAL )
        python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --only MC13TeV_TTJets --exactonly --skipexisting;
        ;;

    FULLSELSYST )
        python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${run_summaryeosdir}_expsyst --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --only MC13TeV_TTJets --systVar all --exactonly;
        ;;

    FULLSELGEN_OLD )
        #ttbar GEN samples
        python scripts/runLocalAnalysis.py -i /eos/user/m/mseidel/ReReco2016/b312177_merged -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --farmappendix samplesGEN;
        ;;

    FULLSELGEN )
        #ttbar GEN samples
        python scripts/runLocalAnalysis.py -i /store/group/cmst3/user/mseidel/ReReco2016/b312177_GEN_merged -q ${queue} -o ${run_summaryeosdir} --era era2016 -m TOPJetShape::RunTopJetShape --skipexisting --farmappendix samplesGEN;
        ;;

    MERGE )
        # options: -T outdir -f
        python scripts/mergeOutputs.py ${summaryeosdir} True - False;
        ;;

    PLOTSEL )
        rm -r plots
        commonOpts="-i ${summaryeosdir} -j data/era2016/samples.json,data/era2016/qcd_samples.json --systJson data/era2016/syst_samples.json,data/era2016/expsyst_samples.json -l ${lumi} --mcUnc ${lumiUnc} --rebin 1"
        python scripts/plotter.py ${commonOpts} --outDir plots;
        ;;

    TESTPLOTSEL )
        commonOpts="-i ${summaryeosdir} -j data/era2016/samples.json,data/era2016/qcd_samples.json --systJson data/era2016/syst_samples.json,data/era2016/expsyst_samples.json -l ${lumi}"
        python scripts/plotter.py ${commonOpts} --outDir plots/test --only L4_1l4j2b2w_njets,L4_1l4j2b2w_nvtx,js_tau32_charged,js_mult_charged,L4_1l4j2b2w_l0eta,L2_1l4j_nbjets,L3_1l4j2b_nwjets;
        #python scripts/plotter.py ${commonOpts} --outDir plots/test --only L4_1l4j2b2w_nvtx;
        ;;

    WWWSEL )
        rm -r ${wwwdir}/sel
        mkdir -p ${wwwdir}/sel
        cp plots/*.{png,pdf} ${wwwdir}/sel
        cp test/index.php ${wwwdir}/sel
        ;;

    BINNING )
        python test/TopJSAnalysis/optimizeUnfoldingMatrix.py -i eos --obs all
        ;;

    WWWBINNING )
        rm -r ${wwwdir}/binning
        mkdir -p ${wwwdir}/binning
        cp unfolding/optimize/*.{png,pdf} ${wwwdir}/binning
        cp test/index.php ${wwwdir}/binning
        ;;

    FILL )
        cd batch;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -q workday -i ${summaryeosdir}/Chunks/ --skip "MC13TeV_TTJets,MC13TeV_TTJets_herwigpp_asfsr*,MC13TeV_TTJets_pythia8_asfsr*,MC13TeV_TTJets_jec_*" --skipexisting --farmappendix fill;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -q workday -i ${summaryeosdir}/Chunks/ --only "MC13TeV_TTJets" --skipexisting --weights pu_down,pu_up,leptrig_up,leptrig_down,lepsel_up,lepsel_down,bfrag_up,bfrag_down,bfrag_peterson,slbr_up,slbr_down,top_pt,id1005muR2muF2hdampmt272.7225,id1009muR0.5muF0.5hdampmt272.7225,id3001PDFset13100,id4001PDFset25200 --farmappendix fill_TTJets;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -q workday -i ${summaryeosdir}/Chunks/ --only "MC13TeV_TTJets_herwigpp_asfsr*,MC13TeV_TTJets_pythia8_asfsr*" --skipexisting --farmappendix fill_asfsr;
        python ../test/TopJSAnalysis/fillUnfoldingMatrix.py -q workday -i ${summaryeosdir}/Chunks/ --only "MC13TeV_TTJets_jec_*" --skipexisting --farmappendix fill_jec;
        cd -;
        ;;

    MERGEFILL )
        # options: -T outdir -f
        ./scripts/mergeOutputs.py unfolding/fill True - False
        ;;

    TOYUNFOLDING )
        mkdir unfolding/toys;
        # mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_02 c1_05 c1_10 c1_20 c2_02 c2_05 c2_10 c2_20 c3_02 c3_05 c3_10 c3_20 m2_b1 n2_b1 n3_b1 m2_b2 n2_b2 n3_b2
        for OBS in mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_00 c1_02 c1_05 c1_10 c1_20 c2_00 c2_02 c2_05 c2_10 c2_20 c3_00 c3_02 c3_05 c3_10 c3_20 m2_b1 n2_b1 n3_b1 m2_b2 n2_b2 n3_b2 nsd
        do
          for RECO in charged all
          do
            for FLAVOR in incl bottom light gluon
            do
              while [ $(jobs | wc -l) -ge 4 ] ; do sleep 1 ; done
              mkdir -p unfolding/toys_farm/${OBS}_${RECO}_${FLAVOR}
              cp test/TopJSAnalysis/testUnfold0Toys.C unfolding/toys_farm/${OBS}_${RECO}_${FLAVOR}
              root -l -b -q "unfolding/toys_farm/${OBS}_${RECO}_${FLAVOR}/testUnfold0Toys.C++(\"${OBS}\", \"${RECO}\", \"${FLAVOR}\", 1000)"&
            done
          done
        done
        ;;

    UNFOLDING )
        for OBS in mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_00 c1_02 c1_05 c1_10 c1_20 c2_00 c2_02 c2_05 c2_10 c2_20 c3_00 c3_02 c3_05 c3_10 c3_20 m2_b1 n2_b1 n3_b1 m2_b2 n2_b2 n3_b2 nsd
        do
          for RECO in charged all
          do
            for FLAVOR in incl bottom light gluon
            do
              while [ $(jobs | wc -l) -ge 4 ] ; do sleep 1 ; done
              python test/TopJSAnalysis/doUnfolding.py --obs ${OBS} --reco ${RECO} --flavor ${FLAVOR} &
            done
          done
        done
        python test/TopJSAnalysis/plotMeanTau.py
        for FLAVOR in incl bottom light gluon
        do
          python test/TopJSAnalysis/plotMeanCvsBeta.py --flavor ${FLAVOR}
        done
        python test/TopJSAnalysis/doCovarianceAndChi2.py
        python test/TopJSAnalysis/doCovarianceAndChi2.py --reco all
        ;;

    WWWUNFOLDING )
        rm -r ${wwwdir}/result
        mkdir -p ${wwwdir}/result
        cp unfolding/result/*.{png,pdf} ${wwwdir}/result
        cp test/index.php ${wwwdir}/result
        ;;

    GENPTPLOTS )
        for FLAVOR in incl bottom light gluon
        do
          python test/TopJSAnalysis/quickPlotMCGenerators.py -n 100000000 --flavor ${FLAVOR} &
        done
        ;;

    FLAVORPLOTS )
        for RECO in charged all
          do
          for OBS in mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_00 c1_02 c1_05 c1_10 c1_20 c2_00 c2_02 c2_05 c2_10 c2_20 c3_00 c3_02 c3_05 c3_10 c3_20 m2_b1 n2_b1 n3_b1 m2_b2 n2_b2 n3_b2 nsd
          do
            python test/TopJSAnalysis/compareFlavors.py --obs ${OBS} --reco ${RECO}
          done
        done
        ;;

    SCANAS )
        for OBS in lowcornew c1all lambda lowcor #c1pert
        do
          for GENERATOR in pythia8 #herwigpp
          do
            for FLAVOR in incl bottom light gluon
            do
              for RECO in charged all
              do
                python test/TopJSAnalysis/scanAlphaS.py --generator ${GENERATOR} --flavor ${FLAVOR} --obs ${OBS} --reco ${RECO}
              done
            done
          done
        done
        ;;
    
    FITAS )
        for OBS in zgdr
        do
          for GENERATOR in pythia8 #herwigpp
          do
            for FLAVOR in incl bottom light gluon
            do
              for RECO in charged all
              do
                python test/TopJSAnalysis/fitAlphaS.py --generator ${GENERATOR} --flavor ${FLAVOR} --obs ${OBS} --reco ${RECO} --cmw
              done
            done
          done
        done
        ;;

    YODA )
        python test/TopJSAnalysis/exportToYoda.py
        ;;

    SENSITIVITY )
        for OBS in mult width ptd ptds ecc tau21 tau32 tau43 zg zgxdr zgdr ga_width ga_lha ga_thrust c1_00 c1_02 c1_05 c1_10 c1_20 c2_00 c2_02 c2_05 c2_10 c2_20 c3_00 c3_02 c3_05 c3_10 c3_20 m2_b1 n2_b1 n3_b1 m2_b2 n2_b2 n3_b2 nsd
        do
          python test/TopJSAnalysis/fitAlphaS.py --generator pythia8 --flavor incl --obs ${OBS}
        done
        ;;

    WWWTUNING )
        rm -r ${wwwdir}/fit
        mkdir -p ${wwwdir}/fit
        cp fit/*.{png,pdf} ${wwwdir}/fit
        cp test/index.php ${wwwdir}/fit
        ;;

    WWW )
        sh scripts/steerTOPJS.sh WWWSEL;
        sh scripts/steerTOPJS.sh WWWBINNING;
        sh scripts/steerTOPJS.sh WWWUNFOLDING;
        sh scripts/steerTOPJS.sh WWWTUNING;
        ;;

esac
