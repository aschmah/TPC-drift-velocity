
// first do
// export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH

R__LOAD_LIBRARY(calib_TPC_vD_t0_C.so);

void run_calib_TPC_vD_t0(TString maininputdir = "/Users/aschmah/alice/TPC_calibration/drift_velocity/download_data/",
                         TString inputfilename = "residuals/alice/data/2022/LHC22m/523309/apass3/0150/o2_ctf_run00523309_orbit0368904604_tf0000000001_epn188/o2tpc_residuals_1660435003162_1660435603148_0_52712.root",
                         Int_t runNumber = 523309)
{

    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1669433390606_1669433990592_632556_685268.root",529663)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1684529925253_1684530525248_843408_1054259.root",536663)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","merge_after_calib_unbinned_529403.root",529403)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/drift_velocity/download_data/residuals/alice/data/2022/LHC22s/529403/apass5/1830/o2_ctf_run00529403_orbit0059839744_tf0000000001_epn205/","o2tpc_residuals_1668791871280_1668792471267_0_52712.root",529403)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1668788039413_1668788639399_0_52712_apass5_with_vD.root",529399)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1666810586060_1666811042049_0_40061_528093.root",528093)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1668510611345_1668511211332_0_52712_529235.root",529235)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1685475045294_1685475645289_421704_632555_537249.root",537249)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/LowIntWithMat/","merge0_10_528675_4kHz.root",528675)
    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/LowIntWithMat/","merge0_10_528675_4kHz_GEO.root",528675)

    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/","o2tpc_residuals_1665629318522_1665629918509_0_52712_527228.root",527228)

    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/drift_velocity/download_data/","o2tpc_residuals_1655115672273_1655116272259_52713_105425_518541.root",518541)

    // .x run_calib_TPC_vD_t0.C("/Users/aschmah/alice/TPC_calibration/drift_velocity/download_data/residuals/","o2tpc_residuals_1665184139340_1665184739327_263565_316277.root",526865)

    //gROOT->ProcessLine(".L calib_TPC_vD_t0.C++");
    printf("run_calib_TPC_vD_t0 started \n");

    gSystem ->Load("calib_TPC_vD_t0_C.so");
    calib_TPC_vD_t0(maininputdir.Data(),inputfilename.Data(),runNumber);
    //gROOT->ProcessLine(".q");
}