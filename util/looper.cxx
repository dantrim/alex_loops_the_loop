// EDM
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "AsgTools/ToolHandle.h"
#include "xAODCore/ShallowCopy.h"

// Base
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "AsgTools/AnaToolHandle.h"

// Muons
#include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"

// Electrons
#include "EgammaAnalysisInterfaces/IEgammaCalibrationAndSmearingTool.h"
#include "EgammaAnalysisInterfaces/IAsgElectronLikelihoodTool.h"

// std/stl
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <stdexcept>
using namespace std;

// ROOT
#include "TFile.h"
#include "TH2F.h"

enum ElectronId {
    LooseLLH=0
    ,MediumLLH
    ,TightLLH
    ,ElectronIdInvalid
};

enum MuonId {
    Tight=0
    ,Medium
    ,Loose
    ,HighPt
    ,LowPt
    ,MuonIdInvalid
};

string ElectronId2Str(ElectronId id_in)
{
    string id = "Invalid";
    switch(id_in) {
        case(ElectronId::LooseLLH) : { id = "LooseBLLHElectron"; break; }
        case(ElectronId::MediumLLH) : { id = "MediumLHElectron"; break; }
        case(ElectronId::TightLLH) : { id = "TightLHElectron"; break; }
        case(ElectronId::ElectronIdInvalid) : { id = "Invalid"; break; }
    }
    return id;
}

string MuonId2Str(MuonId id_in)
{
    string id = "Invalid";
    switch(id_in) {
        case(MuonId::Loose) : { id = "Loose"; break; }
        case(MuonId::Medium): { id = "Medium"; break; }
        case(MuonId::Tight) : { id = "Tight"; break; }
        case(MuonId::LowPt) : { id = "LowPt"; break; }
        case(MuonId::HighPt) : { id = "HighPt"; break; }
        case(MuonId::MuonIdInvalid) : { id = "Invalid"; break; }
    }
    return id;
}

void print_usage(string name)
{
    cout << name << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "   -i|--input          input DAOD file" << endl;
    cout << "   -n|--nevents        number of events to process [default: all]" << endl;
    cout << "   --isatlfast         set if you are running AFII [default: false]" << endl;
    cout << "   -h|--help           print this help message" << endl;
    cout << " Example:" << endl;
    cout << "  " << name << " -i <input-file> -n 500" << endl;
}

TH2F* get_histo(int max, string type)
{
    string histo_titles = "";
    if(type == "el") {
        histo_titles = "Electron WP;;";
    }
    else if(type == "mu") {
        histo_titles = "Muon WP;;";
    }
    string histo_name = "h_" + type;

    TH2F* h = new TH2F(histo_name.c_str(), histo_titles.c_str(), max, 0, max, max, 0, max);
    vector<string> names;
    for(int i = 0; i < max; i++) {
        if(type == "el") {
            names.push_back(ElectronId2Str((ElectronId)i));
        }
        else if(type == "mu") {
            names.push_back(MuonId2Str((MuonId)i));
        }
    }
    cout << "Histo Names: ";
    for(auto name : names) cout << " " << name;
    cout << endl;
    for(int i = 0; i < max; i++) {
        h->GetXaxis()->SetBinLabel(i+1, names.at(i).c_str());
        h->GetYaxis()->SetBinLabel(i+1, names.at(i).c_str());
    }

    return h;
}

int main(int argc, char** argv)
{
    const char * algo_name = argv[0];
    string filename = "";
    int n_to_process = -1;
    bool is_atlfast = false;
    int optin = 1;
    if(argc == 1) {
        print_usage(algo_name);
        return 1;
    }
    while(optin < argc) {
        string in = argv[optin];
        if(in == "-i" || in == "--input") filename = argv[++optin];
        else if(in == "-n" || in == "--nevents") n_to_process = atoi(argv[++optin]);
        else if(in == "--atlfast") is_atlfast = true;
        else if(in == "-h" || in == "--help") { print_usage(algo_name); return 0; }
        else {
            cout << algo_name << " Unknown command line argument '" << in << "' provided, exiting" << endl;
            return 1;
        }
        optin++;
    } // while

    if(filename == "") {
        cout << algo_name << "    ERROR you must provide an input file" << endl;
        return 1;
    }

    RETURN_CHECK(algo_name, xAOD::Init());
    xAOD::TEvent event(xAOD::TEvent::kClassAccess);

    std::unique_ptr<TFile> ifile(TFile::Open(filename.c_str(), "READ"));
    if(!ifile.get() || ifile->IsZombie()) {
        throw std::logic_error("ERROR Could not open input file: " + filename);
    }

    TFile* out_rfile = new TFile("loop_fun.root", "RECREATE");
    TH1F* electron_residual = new TH1F("h_ele_corr", "Electron p_{T} Diff;;", 100, 0, -1);
    electron_residual->GetXaxis()->SetTitle("(p_{T}^{after} - p_{T}^{before}) / p_{T}^{before}");
    TH1F* muon_residual = new TH1F("h_muo_corr", "Muon p_{T} Diff;;", 100, 0, -1);
    muon_residual->GetXaxis()->SetTitle("(p_{T}^{after} - p_{T}^{before}) / p_{T}^{before}");

    TH1F* h_mu_den = new TH1F("h_mu_den", "Muon p_{T}", 20, 0, 100);
    h_mu_den->GetXaxis()->SetTitle("p_{T} [GeV]");
    TH1F* h_mu_num = new TH1F("h_mu_num", "Muon p_{T} Passing LowPt", 20, 0, 100);
    h_mu_num->GetXaxis()->SetTitle("p_{T} [GeV]");

    TH2F* electron_histo = get_histo((int)ElectronId::ElectronIdInvalid, "el");
    TH2F* muon_histo = get_histo((int)MuonId::MuonIdInvalid, "mu");

    asg::AnaToolHandle<CP::IEgammaCalibrationAndSmearingTool> egamma_calib_tool;
    vector<asg::AnaToolHandle<IAsgElectronLikelihoodTool>> electron_id_tools;

    asg::AnaToolHandle<CP::IMuonCalibrationAndSmearingTool> muon_calib_tool;
    vector<asg::AnaToolHandle<CP::IMuonSelectionTool>> muon_id_tools;

    // initialize the tools

    stringstream tool_name;
    stringstream type_name;
    for(int i = 0; i < MuonId::MuonIdInvalid; i++) {
        asg::AnaToolHandle<CP::IMuonSelectionTool> tool("");
        muon_id_tools.push_back(tool);
        tool_name.str("");
        type_name.str("");
        string mu_qual = MuonId2Str((MuonId)i);
        tool_name << "MuonSelectionTool_" + mu_qual; 
        type_name << "CP::MuonSelectionTool/" + tool_name.str();
        cout << algo_name << "   Initializing tool: " << type_name.str() << endl;
        muon_id_tools.at(i).setTypeAndName(type_name.str());
        RETURN_CHECK(algo_name, muon_id_tools.at(i).setProperty("MaxEta", 2.7));
        RETURN_CHECK(algo_name, muon_id_tools.at(i).setProperty("MuQuality", i));
        RETURN_CHECK(algo_name, muon_id_tools.at(i).retrieve());
    } // i

    tool_name.str("");
    type_name.str("");
    for(int i = 0; i < ElectronId::ElectronIdInvalid; i++) {
        asg::AnaToolHandle<IAsgElectronLikelihoodTool> tool("");
        electron_id_tools.push_back(tool);
        tool_name.str("");
        type_name.str("");
        tool_name << "ElectronLLHTool_" + ElectronId2Str((ElectronId)i);
        type_name << "AsgElectronLikelihoodTool/" + tool_name.str();
        cout << algo_name << "   Initializing tool: " << type_name.str() << endl;
        electron_id_tools.at(i).setTypeAndName(type_name.str());
        RETURN_CHECK(algo_name, electron_id_tools.at(i).setProperty("WorkingPoint", ElectronId2Str((ElectronId)i)));
        RETURN_CHECK(algo_name, electron_id_tools.at(i).retrieve());
    } // i

    // calibration tools
    tool_name.str("");
    type_name.str("");
    type_name << "CP::MuonCalibrationPeriodTool/MuonCalibrationAndSmearingTool";
    muon_calib_tool.setTypeAndName(type_name.str());
    RETURN_CHECK(algo_name, muon_calib_tool.retrieve());

    tool_name.str("");
    type_name.str("");
    type_name << "CP::EgammaCalibrationAndSmearingTool/EgammaCalibrationAndSmearingTool";
    egamma_calib_tool.setTypeAndName(type_name.str());
    RETURN_CHECK(algo_name, egamma_calib_tool.setProperty("ESModel", "es2017_R21_v0"));
    RETURN_CHECK(algo_name, egamma_calib_tool.setProperty("decorrelationModel", "1NP_v1"));
    RETURN_CHECK(algo_name, egamma_calib_tool.setProperty("useAFII", is_atlfast));
    RETURN_CHECK(algo_name, egamma_calib_tool.retrieve());

    // start looping
    RETURN_CHECK(algo_name, event.readFrom(ifile.get()));
    const unsigned long long n_entries = event.getEntries();
    cout << algo_name << "    Input file has " << n_entries << " events";
    if(n_to_process < 0) n_to_process = n_entries;
    cout << ", will read " << n_to_process << " events" << endl;

    int n_mu = 0;
    int n_el = 0;

    for(unsigned long long entry = 0; entry < (unsigned)n_to_process; ++entry) {

        bool ok = event.getEntry(entry) >= 0;
        if(!ok) throw std::logic_error("ERROR getEntry failed");

        if(entry % 500 == 0) {
            cout << algo_name << "   Processing entry " << entry << endl;
        }

        // get the object containers
        const xAOD::ElectronContainer* in_electrons = 0;
        const xAOD::MuonContainer* in_muons = 0;

        RETURN_CHECK(algo_name, event.retrieve(in_electrons, "Electrons"));
        RETURN_CHECK(algo_name, event.retrieve(in_muons, "Muons"));

        std::pair<xAOD::ElectronContainer*, xAOD::ShallowAuxContainer*> shallow_ele = xAOD::shallowCopyContainer(*in_electrons);
        std::pair<xAOD::MuonContainer*, xAOD::ShallowAuxContainer*> shallow_muo = xAOD::shallowCopyContainer(*in_muons);
        auto electrons = shallow_ele.first;
        auto muons = shallow_muo.first;

        // loop over electrons
        for(const auto & electron : *electrons) {
            if(electron->pt() < 4e3) continue;
            if(!(fabs(electron->eta()) < 2.47)) continue;
            if(!electron->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)) continue;
            float pt_b4 = electron->pt() * 1e-3;
            if(egamma_calib_tool->applyCorrection(*electron) != CP::CorrectionCode::Ok) continue;
            float pt_after = electron->pt() * 1e-3;
            float diff = (pt_after - pt_b4) / pt_b4;
            electron_residual->Fill(diff);
            n_el++;

            // loop in a loop!
            for(int i = 0; i < (int)ElectronId::ElectronIdInvalid; i++) {
                bool pass_i = electron_id_tools.at(i)->accept(electron);
                if(!pass_i) continue;
                for(int j = 0; j < (int)ElectronId::ElectronIdInvalid; j++) {
                    bool pass_j = electron_id_tools.at(j)->accept(electron);
                    if(pass_j) electron_histo->Fill(i,j);
                } // j
            } // i
        } // electron loop

        for(const auto & muon : *muons) {
            if(muon->pt() < 3e3) continue;
            if(!(fabs(muon->eta()) < 2.5)) continue;
            float pt_b4 = muon->pt() * 1e-3;
            if(!(muon_calib_tool->applyCorrection(*muon) == CP::CorrectionCode::Ok)) continue;
            float pt_after = muon->pt() * 1e-3;
            float diff = (pt_after - pt_b4) / pt_b4;
            muon_residual->Fill(diff);
            n_mu++;

            h_mu_den->Fill(muon->pt() * 1e-3);

            // another loop in a loop!
            for(int i = 0; i < (int)MuonId::MuonIdInvalid; i++) {
                if( i == (int)MuonId::LowPt && muon_id_tools.at(i)->accept(*muon) ) h_mu_num->Fill(muon->pt() * 1e-3);
                bool pass_i = muon_id_tools.at(i)->accept(*muon);
                if(!pass_i) continue;
                for(int j = 0; j < (int)MuonId::MuonIdInvalid; j++) {
                    bool pass_j = muon_id_tools.at(j)->accept(*muon);
                    if(pass_j) muon_histo->Fill(i,j);
                } // j
            } // i
        } // muon loop
    } // entry


    cout << "--------------------------------------" << endl;
    cout << " done looping" << endl;
    cout << algo_name << "    Electron histo has " << n_el << " entries" << endl;
    cout << algo_name << "    Muon histo has " << n_mu << " entries" << endl;


    for(int i = 0; i < (int)MuonId::MuonIdInvalid; i++) {
        for(int j = 0; j < (int)MuonId::MuonIdInvalid; j++) {
            //int bin = muon_histo->GetBin(i,j);
            auto content = muon_histo->GetBinContent(i+1,j+1);
            cout << "Muon " << i << " " << j << " = "  << (content / n_mu) << endl;
            muon_histo->SetBinContent(i+1,j+1, content / n_mu);
        }
    }
    for(int i = 0; i < (int)ElectronId::ElectronIdInvalid; i++) {
        for(int j = 0; j < (int)ElectronId::ElectronIdInvalid; j++) {
            //int bin = electron_histo->GetBin(i,j);
            auto content = electron_histo->GetBinContent(i+1,j+1);
            cout << "Electron " << i << " " << j << " = " << (content / n_el) << endl;
            electron_histo->SetBinContent(i+1,j+1,content / n_el);
        }
    }

    out_rfile->cd();
    electron_histo->Write();
    muon_histo->Write();
    electron_residual->Write();
    muon_residual->Write();
    h_mu_num->Write();
    h_mu_den->Write();

    // clean up memory.. meh

    return 0;
}
