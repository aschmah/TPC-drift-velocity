#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "Framework/Logger.h"
#include "CommonConstants/LHCConstants.h"
#include "SpacePoints/SpacePointsCalibConfParam.h"
#include "SpacePoints/TrackResiduals.h"
#include "SpacePoints/TrackInterpolation.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsTPC/Defs.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsTPC/VDriftCorrFact.h"

#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGrid.h>
#include <TH2.h>
#include <TF1.h>
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom3.h"
#include "TGraph.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <array>
#include <random>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

#include "BC_shifts.h"

#else

#error This macro must run in compiled mode

#endif

using namespace o2::tpc;
using GID = o2::dataformats::GlobalTrackID;

static const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
static const Double_t kAlmost0=Double_t(FLT_MIN);
static const Double_t kB2C=-0.299792458e-3;
static Float_t fHelix[9];
static Float_t track_pos[3];
static char NoP[50];

//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,
                      Float_t size=0.06,Int_t color=1,Float_t angle=0.0,
                      Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color
    // align: 1 left aligned, 32, right aligned

    // align = 10*HorizontalAlign + VerticalAlign
    // For horizontal alignment the following convention applies:
    // 1=left adjusted, 2=centered, 3=right adjusted
    // For vertical alignment the following convention applies:
    // 1=bottom adjusted, 2=centered, 3=top adjusted

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------


void set_helix(Float_t x, Float_t alpha, Float_t param[5], Float_t Bz)
{
    // http://alidoc.cern.ch/AliRoot/master/_ali_helix_8cxx_source.html
    // AliHelix::AliHelix(const AliExternalTrackParam &t)

    Float_t cs,sn;
    for(Int_t i = 0; i < 5; i++)
    {
        fHelix[i] = param[i];
    }

    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...
    fHelix[4]=fHelix[4]/(-1000/0.299792458/Bz);    // C
    //  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
    cs=cosf(alpha); sn=sinf(alpha); // RS use float versions: factor 2 in CPU speed

    Float_t xc, yc, rc;
    rc  =  1/fHelix[4];
    xc  =  x-fHelix[2]*rc;
    Float_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
    if (dummy<0)
    {
        dummy = 0;
    }
    yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];

    fHelix[6] = xc*cs - yc*sn;
    fHelix[7] = xc*sn + yc*cs;
    fHelix[8] =  TMath::Abs(rc);
    //
    //
    fHelix[5]=x*cs - fHelix[0]*sn;            // x0
    fHelix[0]=x*sn + fHelix[0]*cs;            // y0
    //fHelix[1]=                               // z0
    //  fHelix[2]=TMath::ASin(fHelix[2]) + alpha; // phi0
    float sphi = TMath::Abs(fHelix[2])<kAlmost1 ? asinf(fHelix[2]) : TMath::Sign(TMath::Pi()/2, fHelix[2]);
    fHelix[2]=sphi + alpha; // phi0 // RS : use float version
    //fHelix[3]=                               // tgl
    //
    //
    fHelix[5]   = fHelix[6];
    fHelix[0]   = fHelix[7];
}

void evaluate_helix(Float_t t,Float_t r[3])
{
    float phase=fHelix[4]*t+fHelix[2];
    Float_t sn=sinf(phase), cs=cosf(phase);
    //  Float_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

    r[0] = fHelix[5] + sn/fHelix[4];
    r[1] = fHelix[0] - cs/fHelix[4];
    r[2] = fHelix[1] + fHelix[3]*t;
}


TLorentzVector get_TLV_helix(Float_t Bz, Double_t& pt)
{
    TLorentzVector TLV_helix;
    Double_t charge = 1.0;
    Double_t conversion = -1000/0.299792458/Bz;
    pt = charge/(conversion*fHelix[4]);

    //printf("pt: %4.3f \n",pt);
    //fHelix[4] = charge/(conversion*pt); // C

    Float_t track_posA[3];
    evaluate_helix(-85.0,track_posA);
    Float_t track_posB[3];
    evaluate_helix(-85.0+0.1,track_posB);
    TVector3 TV3_dir;
    TV3_dir.SetXYZ(track_posB[0]-track_posA[0],track_posB[1]-track_posA[1],track_posB[2]-track_posA[2]);
    //printf(" \n");
    //printf("posA: {%4.3f, %4.3f}, posB: {%4.3f, %4.3f}, dirA: {%4.3f, %4.3f} \n",track_posA[0],track_posA[1],track_posB[0],track_posB[1],TV3_dir.X(),TV3_dir.Y());
    if(TV3_dir.Mag() > 0.0) TV3_dir *= 1.0/TV3_dir.Mag();
    //printf("dirB: {%4.3f, %4.3f} \n",TV3_dir.X(),TV3_dir.Y());
    Double_t perp = TV3_dir.Perp();
    //printf("perp: %4.3f, pt: %4.3f \n",perp,pt);
    if(perp > 0.0) TV3_dir *= fabs(pt)/perp;
    //printf("dirC: {%4.3f, %4.3f} \n",TV3_dir.X(),TV3_dir.Y());
    //printf("pt: %4.3f, perp: %4.3f \n",pt,TV3_dir.Perp());

    TLV_helix.SetXYZM(TV3_dir.X(),TV3_dir.Y(),TV3_dir.Z(),0.139);
    return TLV_helix;
}


Float_t get_helix_cluster_residuum(Float_t clus_x, Float_t clus_y, Float_t clus_z, Float_t track_path_start, Float_t delta_pos[3], Float_t pos_helix[3], Float_t delta_track_path)
{
    Float_t track_pos[3] = {0.0,0.0,0.0};
    Float_t track_pos_previous[3] = {0.0,0.0,0.0};
    Float_t delta_clus_previous = 99999.0;
    Float_t track_path_minimum = 0.0;
    for(Float_t track_path = track_path_start; fabs(track_path) < 350.0; track_path += delta_track_path)
    {
        evaluate_helix(track_path,track_pos);
        Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);

        //printf("track_path: %4.3f, radius_track: %4.3f \n",track_path,radius_track);

        if(radius_track > 250.0) break;

        Double_t delta_clus = TMath::Sqrt(TMath::Power(clus_x - track_pos[0],2) + TMath::Power(clus_y - track_pos[1],2) +TMath::Power(clus_z - track_pos[2],2));
        if(delta_clus > delta_clus_previous) // minimum reached
        {
            delta_pos[0] = track_pos_previous[0] - clus_x;
            delta_pos[1] = track_pos_previous[1] - clus_y;
            delta_pos[2] = track_pos_previous[2] - clus_z;

            pos_helix[0] = track_pos_previous[0];
            pos_helix[1] = track_pos_previous[1];
            pos_helix[2] = track_pos_previous[2];

            //printf("  --> track_path: %4.3f, delta_clus: %4.3f, delta: {%4.3f, %4.3f, %4.3f} \n",track_path,delta_clus,delta_pos[0],delta_pos[1],delta_pos[2]);

            return track_path_minimum;
        }
        else
        {
            delta_clus_previous = delta_clus;
            track_pos_previous[0] = track_pos[0];
            track_pos_previous[1] = track_pos[1];
            track_pos_previous[2] = track_pos[2];
            track_path_minimum = track_path;
        }
    }

    return track_path_minimum;
}


Float_t get_helix_cluster_residuum_2D(Float_t clus_x, Float_t clus_y, Float_t track_path_start, Float_t delta_pos[3], Float_t pos_helix[3], Float_t delta_track_path, Double_t& delta_clus)
{
    Float_t track_pos[3] = {0.0,0.0,0.0};
    Float_t track_pos_previous[3] = {0.0,0.0,0.0};
    Float_t delta_clus_previous = 99999.0;
    Float_t track_path_minimum = 0.0;
    for(Float_t track_path = track_path_start; fabs(track_path) < 350.0; track_path += delta_track_path)
    {
        evaluate_helix(track_path,track_pos);
        Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);

        //printf("track_path: %4.3f, radius_track: %4.3f \n",track_path,radius_track);

        if(radius_track > 250.0) break;

        delta_clus = TMath::Sqrt(TMath::Power(clus_x - track_pos[0],2) + TMath::Power(clus_y - track_pos[1],2));
        if(delta_clus > delta_clus_previous) // minimum reached
        {
            delta_pos[0] = track_pos_previous[0] - clus_x;
            delta_pos[1] = track_pos_previous[1] - clus_y;

            pos_helix[0] = track_pos_previous[0];
            pos_helix[1] = track_pos_previous[1];
            pos_helix[2] = track_pos_previous[2];

            //printf("  --> track_path: %4.3f, delta_clus: %4.3f, delta: {%4.3f, %4.3f, %4.3f} \n",track_path,delta_clus,delta_pos[0],delta_pos[1],delta_pos[2]);

            return track_path_minimum;
        }
        else
        {
            delta_clus_previous = delta_clus;
            track_pos_previous[0] = track_pos[0];
            track_pos_previous[1] = track_pos[1];
            track_pos_previous[2] = track_pos[2];
            track_path_minimum = track_path;
        }
    }

    return track_path_minimum;
}




void calib_TPC_vD_t0(TString maininputdir = "/Users/aschmah/alice/TPC_calibration/drift_velocity/download_data/", TString inputfilename = "residuals/alice/data/2022/LHC22m/523309/apass3/0150/o2_ctf_run00523309_orbit0368904604_tf0000000001_epn188/o2tpc_residuals_1660435003162_1660435603148_0_52712.root", Int_t runNumber = 523309)
{

    // calib_TPC_vD_t0\(\"/Users/aschmah/alice/TPC_calibration/Data/unbinned/Data/\",\"o2tpc_residuals_1660430573767_1660431173753_0_52712_apass4_v3.root","523308\"\)

    // 523308, 523897
    SetRootGraphicStyle();

    vector<Int_t> vec_BCshift;
    vector<Int_t> vec_runs;
    BC_shifts(vec_BCshift,vec_runs);

    Double_t t0_offset_CCDB_used_in_BC = 0.0;
    std::vector<Int_t>::iterator itterator;
    itterator = std::find(vec_runs.begin(), vec_runs.end(), runNumber);
    if(itterator != vec_runs.end() )
    {
        Int_t index = itterator - vec_runs.begin();
        t0_offset_CCDB_used_in_BC = vec_BCshift[index];
        printf("Run number is in BC list, index: %d, BC shift: %4.1f \n",index,t0_offset_CCDB_used_in_BC);
    }
    else
    {
        printf("Run number is not in BC list \n");
    }


    vector<TH2D*> h2D_dz_vs_z_various;
    vector<TH2D*> h2D_dz_vs_z_various_trunc;
    vector<TProfile*> TP_dz_vs_z_various_trunc;

    TString HistName;
    TF1* func_PolyFitFunc = new TF1("func_PolyFitFunc",PolyFitFunc,-300,300,6);

    h2D_dz_vs_z_various.resize(9); // TRD+ITS all, TRD+ITS low R, TRD+ITS high R, ITS only all, ITS only low R, ITS only high R, ITS only all + TOF, ITS only low R + TOF, ITS only high R + TOF
    h2D_dz_vs_z_various_trunc.resize(9); // TRD+ITS all, TRD+ITS low R, TRD+ITS high R, ITS only all, ITS only low R, ITS only high R, ITS only all + TOF, ITS only low R + TOF, ITS only high R + TOF
    TP_dz_vs_z_various_trunc.resize(9); // TRD+ITS all, TRD+ITS low R, TRD+ITS high R, ITS only all, ITS only low R, ITS only high R, ITS only all + TOF, ITS only low R + TOF, ITS only high R + TOF
    for(Int_t iPad = 0; iPad < 9; iPad++)
    {
        HistName = "h2D_dz_vs_z_various_";
        HistName += iPad;
        h2D_dz_vs_z_various[iPad] = new TH2D(HistName.Data(),HistName.Data(),200,-250,250,200,-20,20);

        HistName = "h2D_dz_vs_z_various_trunc_";
        HistName += iPad;
        h2D_dz_vs_z_various_trunc[iPad] = new TH2D(HistName.Data(),HistName.Data(),200,-250,250,200,-20,20);

        HistName = "TP_dz_vs_z_various_trunc_";
        HistName += iPad;
        TP_dz_vs_z_various_trunc[iPad] = new TProfile(HistName.Data(),HistName.Data(),200,-250,250);
    }



    TString fullinputfilename = maininputdir;
    fullinputfilename += inputfilename;
    TFile* inputFile = TFile::Open(fullinputfilename.Data());

  
    // Get CCDB objects
    auto& ccdbmgr = o2::ccdb::BasicCCDBManager::instance();
    ccdbmgr.setURL("https://alice-ccdb.cern.ch");
    auto runDuration = ccdbmgr.getRunDuration(runNumber);
    auto tRun = runDuration.first + (runDuration.second - runDuration.first) / 2; // time stamp for the middle of the run duration
    printf("time stamp (middle): %lld, start: %lld \n",tRun,runDuration.first);
    ccdbmgr.setTimestamp(tRun);
    ccdbmgr.get<TGeoManager>("GLO/Config/GeometryAligned");
    auto magField = ccdbmgr.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
    o2::base::Propagator::initFieldFromGRP(magField);
    auto propagator = o2::base::Propagator::Instance();
    Float_t B_field = propagator->getNominalBz();
    const o2::base::MatLayerCylSet* matLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdbmgr.get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    propagator->setMatLUT(matLut);
    printf("B_field: %4.5f \n",B_field);

    auto mMatCorrNO  = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    auto mMatCorrON  = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    auto mMatCorrGEO = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    o2::track::TrackParametrization<float>::dim3_t xyz_glo_pos;
    o2::track::TrackParCov trkWork;

    // CTP orbit reset time
    auto orbitResetTimeNS = ccdbmgr.get<std::vector<int64_t>>("CTP/Calib/OrbitReset");
    int64_t orbitResetTimeMS = (*orbitResetTimeNS)[0] * 1e-3;
    LOGP(info, "Orbit reset time in MS is {}", orbitResetTimeMS);

    std::unique_ptr<TTree> treeUnbinnedResiduals;
    std::unique_ptr<TTree> treeTrackData;
    std::unique_ptr<TTree> treeRecords;
    treeUnbinnedResiduals.reset(nullptr);
    treeUnbinnedResiduals.reset((TTree*)inputFile->Get("unbinnedResid"));
    std::vector<UnbinnedResid> unbinnedResiduals, *unbinnedResidualsPtr = &unbinnedResiduals; // unbinned residuals input
    std::vector<TrackDataCompact> trackRefs, *trackRefsPtr = &trackRefs;                      // the track references for unbinned residuals
    std::vector<TrackData> trackData, *trackDataPtr = &trackData;                             // additional track information (chi2, nClusters, track parameters)
    std::vector<uint32_t> orbits, *orbitsPtr = &orbits;                                       // first orbit for each TF in the input data
    treeUnbinnedResiduals->SetBranchAddress("res", &unbinnedResidualsPtr);
    treeUnbinnedResiduals->SetBranchAddress("trackInfo", &trackRefsPtr);
    treeTrackData.reset((TTree*)inputFile->Get("trackData"));
    treeTrackData->SetBranchAddress("trk", &trackDataPtr);

    Long64_t nentries_trackdata = treeTrackData->GetEntries();
    printf("nentries_trackdata: %lld \n",nentries_trackdata);

    const SpacePointsCalibConfParam& params = SpacePointsCalibConfParam::Instance();


    int64_t startTime = -1;
    int64_t stopTime  = -1;
    treeRecords.reset((TTree*)inputFile->Get("records"));
    treeRecords->SetBranchAddress("firstTForbit", &orbitsPtr);
    treeRecords->GetEntry(0); // per input file there is only a single entry in the tree
    Int_t nentries_records = treeRecords->GetEntries();
    printf("nentries_records: %d \n",nentries_records);
    if(startTime < 0)
    {
        uint32_t minFirstOrbit = -1;
        uint32_t minLastOrbit  = -1;
        Int_t n_orbits = 0;
        for(auto orbit : orbits)
        {
            printf("orbit: %d \n",orbit);
            if(orbit < minFirstOrbit)
            {
                minFirstOrbit = orbit;
            }
            if(orbit > minLastOrbit || n_orbits == 0)
            {
                minLastOrbit = orbit;
            }
            n_orbits++;
        }
        startTime = orbitResetTimeMS + minFirstOrbit * o2::constants::lhc::LHCOrbitMUS * 1.e-3;
        stopTime  = orbitResetTimeMS + minLastOrbit  * o2::constants::lhc::LHCOrbitMUS * 1.e-3;
        printf("startTime: %lld, minFirstOrbit: %d, stopTime: %lld, minLastOrbit: %d \n",startTime,minFirstOrbit,stopTime,minLastOrbit);
    }

    auto vDriftTgl = ccdbmgr.getForTimeStamp<o2::tpc::VDriftCorrFact>("TPC/Calib/VDriftTgl", startTime);
    Double_t refVDrift_CCDB  = (Double_t)vDriftTgl->refVDrift;
    Double_t vD_CCDB_used    = (Double_t)vDriftTgl->getVDrift();
    Double_t t0_offset_CCDB_used_in_mus = (Double_t)vDriftTgl->getTimeOffset();
    printf("refVDrift_CCDB: %4.5f, vD_CCDB_used: %4.5f, t0_offset_CCDB_used_in_mus: %4.5f  \n",refVDrift_CCDB,vD_CCDB_used,t0_offset_CCDB_used_in_mus);

    Int_t entries_TRD = 0;
    Int_t entries_TOF = 0;


    for(int iEntry = 0; iEntry < treeUnbinnedResiduals->GetEntries(); ++iEntry)
    {
        treeUnbinnedResiduals->GetEntry(iEntry);
        treeTrackData        ->GetEntry(iEntry);

        for (const auto& res : unbinnedResiduals)
        {
            Float_t dy = res.dy * param::MaxResid / 0x7fff;
            Float_t dz = res.dz * param::MaxResid / 0x7fff;
            //printf("entry: %d, sector: %d, dy: %4.3f, dz: %4.3f \n",iEntry,(Int_t)res.sec,dy,dz);
            //hSector->Fill(res.sec);
        }



        auto nTracks = trackRefs.size();
        for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
        {
            const auto& trkInfo = trackRefs[iTrack];
            const auto& trk     = trackData[iTrack];
            //printf("iTrack: %d \n",(Int_t)iTrack);
            Int_t N_cls_ITS     = trk.nClsITS;
            UShort_t N_TrkltsTRD   = trk.nTrkltsTRD;
            UShort_t TOF_available = trk.clAvailTOF;
            Float_t  chi2TPC       = trk.chi2TPC; // peak at 320, max ~ 500

            // ReconstructionDataFormats/TrackParametrization.h

            Float_t  track_x       = trk.par.getX();
            Float_t  track_alpha   = trk.par.getAlpha();
            Float_t  track_p[5];
            for(Int_t i_p = 0; i_p < 5; i_p++)
            {
                track_p[i_p]   = trk.par.getParam(i_p);
                trkWork.setParam(track_p[i_p],i_p);
            }

            trkWork.setAlpha(track_alpha);
            trkWork.setX(track_x);


            set_helix(track_x,track_alpha,track_p,B_field);

            Double_t pT_track_charge = 0.0;
            TLorentzVector TLV_helix_prim = get_TLV_helix(B_field,pT_track_charge);
            Double_t pt_track  = TLV_helix_prim.Pt();
            Double_t eta_track = TLV_helix_prim.Eta();
            Double_t phi_track = TLV_helix_prim.Phi();
            Double_t px_track  = TLV_helix_prim.Px();
            Double_t py_track  = TLV_helix_prim.Py();
            Double_t pz_track  = TLV_helix_prim.Pz();

            Float_t delta_pos_beam[3];
            Float_t pos_beam[3];
            //Float_t track_path_beam_axis = get_helix_cluster_residuum(0.0,0.0,0.0,0.0,delta_pos_beam,pos_beam,-0.1);
            //evaluate_helix(track_path_beam_axis,pos_beam);
            //Float_t dca_2D = TMath::Sqrt(pos_beam[0]*pos_beam[0] + pos_beam[1]*pos_beam[1]);

            propagator->PropagateToXBxByBz(trkWork, 0.0, 0.95, 2.0, o2::base::Propagator::MatCorrType::USEMatCorrTGeo); // USEMatCorrTGeo, USEMatCorrLUT, USEMatCorrNONE
            trkWork.getXYZGlo(xyz_glo_pos);

            Float_t dca_2D = TMath::Sqrt(xyz_glo_pos[0]*xyz_glo_pos[0] + xyz_glo_pos[1]*xyz_glo_pos[1]);

            //if(TOF_available) printf("N_cls_ITS: %d, TOF_available: %d, N_TrkltsTRD: %d, dca_2D: %4.3f, eta: %4.3f, pt: %4.3f  \n",N_cls_ITS,TOF_available,N_TrkltsTRD,dca_2D,eta_track,pt_track);

            for(unsigned int i = trkInfo.idxFirstResidual; i < trkInfo.idxFirstResidual + trkInfo.nResiduals; ++i)
            {
                const auto& residIn = unbinnedResiduals[i];
                float xPos = param::RowX[residIn.row];
                float radius = xPos;
                float yPos = residIn.y * param::MaxY / 0x7fff + residIn.dy * param::MaxResid / 0x7fff;
                float zPos = residIn.z * param::MaxZ / 0x7fff + residIn.dz * param::MaxResid / 0x7fff;
                Float_t dy = residIn.dy * param::MaxResid / 0x7fff;
                Float_t dz = residIn.dz * param::MaxResid / 0x7fff;
                //printf("iEntry: %d, track: %d, dy: %4.3f, dz: %4.3f, N_cls_ITS: %d \n",iEntry,(Int_t)iTrack,dy,dz,N_cls_ITS);


                //if(TOF_available) printf("  -> N_cls_ITS: %d, TOF_available: %d, N_TrkltsTRD: %d, dca_2D: %4.3f, eta: %4.3f, pt: %4.3f, radius: %4.3f  \n",N_cls_ITS,TOF_available,N_TrkltsTRD,dca_2D,eta_track,pt_track,radius);


                if(dca_2D < 1.0 && fabs(eta_track) < 0.84 && fabs(pt_track) > 0.7)
                {
                    // with TRD
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 3 && radius > 100.0 && radius < 160.0)
                    //if(N_cls_ITS >= 5 && N_TrkltsTRD >= 3 && radius > 0.0 && radius < 250.0)
                    {
                        h2D_dz_vs_z_various[0] ->Fill(zPos,dz);
                        entries_TRD++;
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 3 && radius < 110.0)
                    {
                        h2D_dz_vs_z_various[1] ->Fill(zPos,dz);
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 3 && radius > 190.0)
                    {
                        h2D_dz_vs_z_various[2] ->Fill(zPos,dz);
                    }


                    // without TRD and without TOF
                    if(N_cls_ITS >= 5 && N_TrkltsTRD <= 0  && radius > 100.0 && radius < 160.0 && TOF_available == 0)
                    //if(N_cls_ITS >= 5 && N_TrkltsTRD <= 0  && radius > 0.0 && radius < 250.0 && TOF_available == 0)
                    {
                        h2D_dz_vs_z_various[3] ->Fill(zPos,dz);
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD <= 0 && TOF_available == 0)
                    {
                        if(radius < 110.0)  h2D_dz_vs_z_various[4] ->Fill(zPos,dz);
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD <= 0 && radius > 190.0 && TOF_available == 0)
                    {
                        h2D_dz_vs_z_various[5] ->Fill(zPos,dz);
                    }


                    // with TOF
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 0 && radius > 100.0 && radius < 160.0 && TOF_available == 1)
                    //if(N_cls_ITS >= 5 && N_TrkltsTRD >= 0 && radius > 0.0 && radius < 250.0 && TOF_available == 1)
                    {
                        h2D_dz_vs_z_various[6] ->Fill(zPos,dz);
                        entries_TOF++;
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 0 && TOF_available == 1)
                    {
                        if(radius < 110.0)  h2D_dz_vs_z_various[7] ->Fill(zPos,dz);
                    }
                    if(N_cls_ITS >= 5 && N_TrkltsTRD >= 0 && radius > 190.0 && TOF_available == 1)
                    {
                        h2D_dz_vs_z_various[8] ->Fill(zPos,dz);
                    }
                }


            }

        }

    }



    //------------------------------------------------------------------------------------------------------------
    TCanvas* can_dz_vs_z_various = new TCanvas("can_dz_vs_z_various","can_dz_vs_z_various",10,10,1450,920);
    can_dz_vs_z_various ->Divide(3,3);
    TString TS_label_dz_vs_z_various[9] = {"ITS+TRD, R < 250 cm","ITS+TRD, R < 110 cm","ITS+TRD, R > 190 cm","ITS w/o TOF, R < 250 cm","ITS w/o TOF, R < 110 cm","ITS w/o TOF, R > 190 cm","ITS w TOF, R < 250 cm","ITS w TOF, R < 110 cm","ITS w TOF, R > 190 cm"};
    for(Int_t iPad = 0; iPad < 9; iPad++)
    {
        can_dz_vs_z_various ->cd(iPad+1);
        can_dz_vs_z_various ->cd(iPad+1)->SetFillColor(10);
        can_dz_vs_z_various ->cd(iPad+1)->SetTopMargin(0.12);
        can_dz_vs_z_various ->cd(iPad+1)->SetBottomMargin(0.15);
        can_dz_vs_z_various ->cd(iPad+1)->SetRightMargin(0.02);
        can_dz_vs_z_various ->cd(iPad+1)->SetLeftMargin(0.15);
        can_dz_vs_z_various ->cd(iPad+1)->SetTicks(1,1);
        can_dz_vs_z_various ->cd(iPad+1)->SetGrid(0,0);
        can_dz_vs_z_various ->cd(iPad+1)->SetLogz(1);

        h2D_dz_vs_z_various[iPad] ->GetXaxis()->CenterTitle();
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->CenterTitle();
        h2D_dz_vs_z_various[iPad] ->SetStats(0);
        h2D_dz_vs_z_various[iPad] ->SetTitle("");
        h2D_dz_vs_z_various[iPad] ->GetXaxis()->SetTitleOffset(1.2);
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetTitleOffset(1.3);
        h2D_dz_vs_z_various[iPad] ->GetXaxis()->SetLabelSize(0.06);
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetLabelSize(0.06);
        h2D_dz_vs_z_various[iPad] ->GetXaxis()->SetTitleSize(0.06);
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetTitleSize(0.06);
        h2D_dz_vs_z_various[iPad] ->GetXaxis()->SetNdivisions(505,'N');
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetNdivisions(505,'N');
        h2D_dz_vs_z_various[iPad] ->GetXaxis()->SetTitle("z (cm)");
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetTitle("dz (cm)");
        h2D_dz_vs_z_various[iPad] ->GetYaxis()->SetRangeUser(-7.0,7.0);
        h2D_dz_vs_z_various[iPad] ->DrawCopy("colz");

        plotTopLegend((char*)TS_label_dz_vs_z_various[iPad].Data(),0.24,0.92,0.055,1,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        // Truncation
        TH1D* h_dummy_proj;
        for(Int_t i_bin_x = 1; i_bin_x <= h2D_dz_vs_z_various[iPad]->GetNbinsX(); i_bin_x++)
        {
            Double_t x_val = h2D_dz_vs_z_various[iPad]->GetXaxis()->GetBinCenter(i_bin_x);
            h_dummy_proj = (TH1D*)h2D_dz_vs_z_various[iPad] ->ProjectionY("dummy",i_bin_x,i_bin_x);
            Double_t rms_proj  = h_dummy_proj ->GetRMS();
            Double_t mean_proj = h_dummy_proj ->GetMean();
            for(Int_t i_bin_y = 1; i_bin_y <= h2D_dz_vs_z_various[iPad]->GetNbinsY(); i_bin_y++)
            {
                Double_t y_val = h2D_dz_vs_z_various[iPad]->GetYaxis()->GetBinCenter(i_bin_y);
                if(fabs(y_val - mean_proj) < 2.0*rms_proj)
                {
                    Double_t bin_cont = h2D_dz_vs_z_various[iPad] ->GetBinContent(i_bin_x,i_bin_y);
                    h2D_dz_vs_z_various_trunc[iPad] ->SetBinContent(i_bin_x,i_bin_y,bin_cont);
                }
            }
            delete h_dummy_proj;
        }



        HistName = "TP_dz_vs_z_various_trunc_";
        HistName += iPad;
        TP_dz_vs_z_various_trunc[iPad] = h2D_dz_vs_z_various_trunc[iPad] ->ProfileX(HistName.Data());



    }


    Double_t z_range_min[2] = {-130.0,25.0};
    Double_t z_range_max[2] = {-25.0,130.0};
    //Double_t z_range_min[2] = {-200.0,20.0};
    //Double_t z_range_max[2] = {-20.0,200.0};
    Int_t z_color_fit[2] = {kRed,kBlue};
    Double_t fit_params[9][2][3] = {0.0}; // iPad, z, par0/par1/chi2

    TCanvas* can_dz_vs_z_various_trunc = new TCanvas("can_dz_vs_z_various_trunc","can_dz_vs_z_various_trunc",10,10,1450,920);
    can_dz_vs_z_various_trunc ->Divide(3,3);
    TString TS_label_dz_vs_z_various_trunc[9] = {"ITS+TRD, 100 < R < 160 cm","ITS+TRD, R < 110 cm","ITS+TRD, R > 190 cm","ITS w/o TOF, 100 < R < 160 cm","ITS w/o TOF, R < 110 cm","ITS w/o TOF, R > 190 cm","ITS w TOF, 100 < R < 160 cm","ITS w TOF, R < 110 cm","ITS w TOF, R > 190 cm"};


    Double_t offset_dz_vs_z = 0.0;
    Double_t slope_dz_vs_z  = 0.0;
    Double_t arr_average_slope[9];
    Double_t arr_average_offset[9];
    for(Int_t iPad = 0; iPad < 9; iPad++)
    {
        can_dz_vs_z_various_trunc ->cd(iPad+1);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetFillColor(10);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetTopMargin(0.12);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetBottomMargin(0.15);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetRightMargin(0.02);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetLeftMargin(0.15);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetTicks(1,1);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetGrid(0,0);
        can_dz_vs_z_various_trunc ->cd(iPad+1)->SetLogz(1);

        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->CenterTitle();
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->CenterTitle();
        h2D_dz_vs_z_various_trunc[iPad] ->SetStats(0);
        h2D_dz_vs_z_various_trunc[iPad] ->SetTitle("");
        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->SetTitleOffset(1.2);
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetTitleOffset(1.3);
        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->SetLabelSize(0.06);
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetLabelSize(0.06);
        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->SetTitleSize(0.06);
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetTitleSize(0.06);
        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->SetNdivisions(505,'N');
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetNdivisions(505,'N');
        h2D_dz_vs_z_various_trunc[iPad] ->GetXaxis()->SetTitle("z (cm)");
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetTitle("dz (cm)");
        h2D_dz_vs_z_various_trunc[iPad] ->GetYaxis()->SetRangeUser(-7.0,7.0);
        h2D_dz_vs_z_various_trunc[iPad] ->DrawCopy("colz");

        plotTopLegend((char*)TS_label_dz_vs_z_various_trunc[iPad].Data(),0.24,0.92,0.055,1,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        TP_dz_vs_z_various_trunc[iPad] ->SetLineColor(kBlack);
        TP_dz_vs_z_various_trunc[iPad] ->SetLineWidth(2);
        TP_dz_vs_z_various_trunc[iPad] ->DrawCopy("hist same");

        for(Int_t i_z = 0; i_z < 2; i_z++)
        {
            for(Int_t i = 0; i < 6; i++)
            {
                func_PolyFitFunc ->SetParameter(i,0.0);
                func_PolyFitFunc ->SetParError(i,0.0);
            }
            func_PolyFitFunc ->SetParameter(0,0.0);
            func_PolyFitFunc ->SetParameter(1,0.0);
            func_PolyFitFunc ->FixParameter(2,0.0);
            func_PolyFitFunc ->FixParameter(3,0.0);
            func_PolyFitFunc ->FixParameter(4,0.0);
            func_PolyFitFunc ->FixParameter(5,0.0);

            TP_dz_vs_z_various_trunc[iPad] ->Fit("func_PolyFitFunc","QWMN","",z_range_min[i_z],z_range_max[i_z]);

            func_PolyFitFunc ->SetLineColor(z_color_fit[i_z]);
            func_PolyFitFunc ->SetLineStyle(1);
            func_PolyFitFunc ->SetLineWidth(4);
            func_PolyFitFunc ->SetRange(z_range_min[i_z],z_range_max[i_z]);

            func_PolyFitFunc ->DrawCopy("same");

            HistName = "p0 = ";
            sprintf(NoP,"%4.5f",(Double_t)func_PolyFitFunc->GetParameter(0));
            HistName += NoP;
            HistName += ", p1 = ";
            sprintf(NoP,"%4.5f",(Double_t)func_PolyFitFunc->GetParameter(1));
            HistName += NoP;
            if(i_z == 0) plotTopLegend((char*)HistName.Data(),0.18,0.3,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(i_z == 1) plotTopLegend((char*)HistName.Data(),0.54,0.7,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

            fit_params[iPad][i_z][0] = (Double_t)func_PolyFitFunc->GetParameter(0);
            fit_params[iPad][i_z][1] = (Double_t)func_PolyFitFunc->GetParameter(1);
            fit_params[iPad][i_z][2] = (Double_t)func_PolyFitFunc->GetChisquare();
        }


        Double_t average_p0 = (-fit_params[iPad][0][0] + fit_params[iPad][1][0])/2.0;
        Double_t average_p1 = (fit_params[iPad][0][1] + fit_params[iPad][1][1])/2.0;
        HistName = "<p0> = ";
        sprintf(NoP,"%4.5f",(Double_t)average_p0);
        HistName += NoP;
        HistName += ", <p1> = ";
        sprintf(NoP,"%4.5f",(Double_t)average_p1);
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.34,0.8,0.045,kGreen+2,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        arr_average_slope[iPad]  = average_p1;
        arr_average_offset[iPad] = average_p0;

        //if(iPad == 6)
        //{
        //    offset_dz_vs_z = average_p0;
        //    slope_dz_vs_z  = average_p1;
        //}

    }

    can_dz_vs_z_various_trunc->cd(1);
    HistName = "run: ";
    HistName += runNumber;
    plotTopLegend((char*)HistName.Data(),0.75,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    can_dz_vs_z_various_trunc->cd(2);
    HistName = "start: ";
    HistName += startTime;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    can_dz_vs_z_various_trunc->cd(3);
    HistName = "stop: ";
    HistName += stopTime;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    //offset_dz_vs_z = 0.5*(arr_average_offset[0] + arr_average_offset[6]);
    //slope_dz_vs_z  = 0.5*(arr_average_slope[0]  + arr_average_slope[6]);

    if(entries_TRD + entries_TOF > 0)
    {
        offset_dz_vs_z = (entries_TRD*arr_average_offset[0] + entries_TOF*arr_average_offset[6])/(entries_TRD + entries_TOF);
        slope_dz_vs_z  = (entries_TRD*arr_average_slope[0]  + entries_TOF*arr_average_slope[6])/(entries_TRD + entries_TOF);
    }

    can_dz_vs_z_various_trunc->cd(4);
    HistName = "<offset>: ";
    sprintf(NoP,"%4.5f",(Double_t)offset_dz_vs_z);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    can_dz_vs_z_various_trunc->cd(5);
    HistName = "<slope>: ";
    sprintf(NoP,"%4.5f",(Double_t)slope_dz_vs_z);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //------------------------------------------------------------------------------------------------------------




    //------------------------------------------------------------------------------------------------------------
    // BC to time in mus:

#if 1
    // values for 523308
    if(runNumber == 523308) // CCDB values were already changed but data is coming from apass3 with old CCDB entries
    {
        t0_offset_CCDB_used_in_BC  = 8.0;
        t0_offset_CCDB_used_in_mus = 0.0;
        vD_CCDB_used               = 2.61894;
    }
    if(
       runNumber == 529397
       || runNumber == 529399
       || runNumber == 529403
       || runNumber == 529414
       || runNumber == 529418
      ) // CCDB values were already changed but data is coming from apass4 with old CCDB entries
    {
        t0_offset_CCDB_used_in_BC  = 86; // was used for apass5 + some vD and t0 from p+p calibrated data
        //t0_offset_CCDB_used_in_BC  = 0;
        t0_offset_CCDB_used_in_mus = 2.4504173;
        vD_CCDB_used               = 2.6891341;
    }
#endif

    //t0_offset_CCDB_used_in_BC  = 0.0; // for apass5 data


    Double_t t0_offset_CCDB_used = -t0_offset_CCDB_used_in_BC * o2::constants::lhc::LHCBunchSpacingNS * 1e-3 - t0_offset_CCDB_used_in_mus;
    printf("t0 offset BC: %4.1f, t0_offset_CCDB_used: %4.5f mus \n",t0_offset_CCDB_used_in_BC,t0_offset_CCDB_used);

    Double_t vD_corr        = vD_CCDB_used/(slope_dz_vs_z + 1.0);
    Double_t delta_t        = (offset_dz_vs_z + 250.0*((vD_CCDB_used/vD_corr) - 1.0))/vD_CCDB_used;
    printf("delta_t: %4.5f mus, t0_offset_CCDB_used_in_BC: %4.5f, t0_offset_CCDB_used_in_mus: %4.5f \n",delta_t,t0_offset_CCDB_used_in_BC,t0_offset_CCDB_used_in_mus);
    //Double_t t0_offset_corr = -(delta_t - t0_offset_CCDB_used); // ???
    Double_t t0_offset_corr = -(delta_t + t0_offset_CCDB_used); // ???

    printf("vD_corr: %4.5f cm/mus \n",vD_corr);
    printf("t0: %4.5f mus \n",t0_offset_corr);

#if 1
    TFile *_file0 = TFile::Open("./snapshot_files/snapshot_dummy.root");
    auto foo = (o2::tpc::VDriftCorrFact*) _file0->Get("ccdb_object");
    static char NoP[50];
    TString outname = "./snapshot_files/snapshot_";
    outname += startTime;
    outname += "_";
    outname += stopTime;
    outname += "_vD_";
    sprintf(NoP,"%4.5f",(Double_t)vD_corr);
    outname += NoP;
    outname += "_t0_";
    sprintf(NoP,"%4.5f",(Double_t)t0_offset_corr);
    outname += NoP;
    outname += "_run_";
    outname += runNumber;
    outname += ".root";
    auto fOut = new TFile(outname, "RECREATE");
    o2::tpc::VDriftCorrFact vdrift_corr = *foo;
    Double_t vD_ref  = vdrift_corr.refVDrift;
    Double_t vD_file = vdrift_corr.getVDrift();
    //vdrift_corr.corrFactErr;

    Double_t corr_factor = vD_corr/vD_ref;
    printf("corr_factor: %4.5f \n",corr_factor);
    vdrift_corr.corrFact = corr_factor;
    Double_t vD_check = vdrift_corr.getVDrift();
    vdrift_corr.timeOffsetCorr = t0_offset_corr;
    printf("vD written to file: %4.5f \n",vD_check);
    fOut->WriteObjectAny(&vdrift_corr, "o2::tpc::VDriftCorrFact","ccdb_object");
    fOut->ls();
    fOut->Close();
#endif

    can_dz_vs_z_various_trunc->cd(7);
    HistName = "vD: ";
    sprintf(NoP,"%4.5f",(Double_t)vD_corr);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    can_dz_vs_z_various_trunc->cd(8);
    HistName = "t0: ";
    sprintf(NoP,"%4.5f",(Double_t)t0_offset_corr);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.7,0.92,0.055,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //------------------------------------------------------------------------------------------------------------



    //------------------------------------------------------------------------------------------------------------
    TString outname_QA_tree = "./QA_tree/LHC22o_vD_t0_";
    outname_QA_tree += startTime;
    outname_QA_tree += "_";
    outname_QA_tree += stopTime;
    outname_QA_tree += "_vD_";
    sprintf(NoP,"%4.5f",(Double_t)vD_corr);
    outname_QA_tree += NoP;
    outname_QA_tree += "_t0_";
    sprintf(NoP,"%4.5f",(Double_t)t0_offset_corr);
    outname_QA_tree += NoP;
    outname_QA_tree += "_run_";
    outname_QA_tree += runNumber;
    outname_QA_tree += "_v2.root";
    TFile* outputfile_vD = new TFile(outname_QA_tree.Data(),"RECREATE");

    Double_t startTime_tree, stopTime_tree, runNumber_tree, t0_offset_CCDB_used_in_BC_tree, vD_CCDB_used_tree, t0_offset_CCDB_used_in_mus_tree;
    Double_t entries_TRD_tree, entries_TOF_tree, chi2TRDnegZ, chi2TRDposZ, chi2TOFnegZ, chi2TOFposZ;

    TTree output_tree("vD_t0_tree","tree vor ALICE TPC vD and t0 calibration");
    output_tree.Branch("vD",&vD_corr,"vD/D");
    output_tree.Branch("t0",&t0_offset_corr,"t0/D");
    output_tree.Branch("tstart",&startTime_tree,"tstart/D");
    output_tree.Branch("tstop",&stopTime_tree,"tstop/D");
    output_tree.Branch("run",&runNumber_tree,"run/D");
    output_tree.Branch("bc",&t0_offset_CCDB_used_in_BC_tree,"bc/D");
    output_tree.Branch("vD_CCDB",&vD_CCDB_used_tree,"vD_CCDB/D");
    output_tree.Branch("t0_CCDB",&t0_offset_CCDB_used_in_mus_tree,"t0_CCDB/D");
    output_tree.Branch("nTRD",&entries_TRD_tree,"nTRD/D");
    output_tree.Branch("nTOF",&entries_TOF_tree,"nTOF/D");
    output_tree.Branch("chi2TRDnegZ",&chi2TRDnegZ,"chi2TRDnegZ/D");
    output_tree.Branch("chi2TRDposZ",&chi2TRDposZ,"chi2TRDposZ/D");
    output_tree.Branch("chi2TOFnegZ",&chi2TOFnegZ,"chi2TOFnegZ/D");
    output_tree.Branch("chi2TOFposZ",&chi2TOFposZ,"chi2TOFposZ/D");


    startTime_tree                  = (Double_t)startTime;
    stopTime_tree                   = (Double_t)stopTime;
    runNumber_tree                  = (Double_t)runNumber;
    t0_offset_CCDB_used_in_BC_tree  = (Double_t)t0_offset_CCDB_used_in_BC;
    vD_CCDB_used_tree               = (Double_t)vD_CCDB_used;
    t0_offset_CCDB_used_in_mus_tree = (Double_t)t0_offset_CCDB_used_in_mus;
    entries_TRD_tree                = (Double_t)entries_TRD;
    entries_TOF_tree                = (Double_t)entries_TOF;
    chi2TRDnegZ                     = (Double_t)fit_params[0][0][2];
    chi2TRDposZ                     = (Double_t)fit_params[0][1][2];
    chi2TOFnegZ                     = (Double_t)fit_params[6][0][2];
    chi2TOFposZ                     = (Double_t)fit_params[6][1][2];

    output_tree.Fill();

    outputfile_vD ->cd();
    output_tree.Write();
    //------------------------------------------------------------------------------------------------------------



    //------------------------------------------------------------------------------------------------------------
    TString outname_QA = "./QA/AC_dz_vs_z_";
    outname_QA += startTime;
    outname_QA += "_";
    outname_QA += stopTime;
    outname_QA += "_vD_";
    sprintf(NoP,"%4.5f",(Double_t)vD_corr);
    outname_QA += NoP;
    outname_QA += "_t0_";
    sprintf(NoP,"%4.5f",(Double_t)t0_offset_corr);
    outname_QA += NoP;
    outname_QA += "_run_";
    outname_QA += runNumber;
    outname_QA += "_v2.png";
    can_dz_vs_z_various_trunc ->SaveAs(outname_QA.Data());
    //------------------------------------------------------------------------------------------------------------

    outputfile_vD ->Close();




}
