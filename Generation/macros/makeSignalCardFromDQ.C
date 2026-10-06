/// Example signal-card writer for Ditto from DQ output

#include "DittoSignal.h"
#include "DittoTune.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

void makeSignalCardFromDQ(const char* tuneFile = "Tuning/tunes/Ditto_tune_pythia8_inel_136tev.root",
                          const char* outputFile = "signal_injected_jpsi_lc.root")
{

  TFile ptjpsi("/tmp/ptjpsi.root");
  //   ptjpsi.ls();
  TH1F* ptjpsi_hist = (TH1F*)ptjpsi.Get("PtMC");
  ptjpsi_hist->SetName("ptjpsi");
  TFile etajpsi("/tmp/etajpsi.root");
  //   etajpsi.ls();
  TH1F* etajpsi_hist = (TH1F*)etajpsi.Get("EtaMC");
  etajpsi_hist->SetName("etajpsi");
  TFile ptlc("/tmp/ptlc.root");
  //   ptlc.ls();
  TH1F* ptlc_hist = (TH1F*)ptlc.Get("PtMC");
  ptlc_hist->SetName("ptlc");
  TFile etalc("/tmp/etalc.root");
  //   etalc.ls();
  TH1F* etalc_hist = (TH1F*)etalc.Get("EtaMC");
  etalc_hist->SetName("etalc");

  auto tune = Ditto::Tune::load(tuneFile);
  if (!tune || tune->mActivityEdges.size() < 2) {
    throw std::runtime_error("Could not load valid Ditto activity edges");
  }

  const auto& activityEdges = tune->mActivityEdges;
  const int nActivity = static_cast<int>(activityEdges.size()) - 1;

  const std::vector<int> pdgs = {443, 4122}; // J/psi, Lambda_c

  const std::vector<double> masses = {3.0969, 2.28646}; // GeV/c^2
  const std::vector<TH1F*> pt_hists = {ptjpsi_hist, ptlc_hist};
  const std::vector<TH1F*> eta_hists = {etajpsi_hist, etalc_hist};

  if (pdgs.size() != masses.size()) {
    throw std::runtime_error("PDG and mass vectors have different sizes");
  }

  TFile* output(TFile::Open(outputFile, "RECREATE"));
  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("Could not create ") + outputFile);
  }

  TH1D* signalMass = Ditto::SignalInjector::makeSignalMassHistogram(*tune, pdgs, masses);
  signalMass->Write();

  TParameter<int> signalDecay("signalDecay", 1);
  signalDecay.Write();

  for (std::size_t iSpecies = 0; iSpecies < pdgs.size(); ++iSpecies) {
    const int pdg = pdgs[iSpecies];
    const double mass = masses[iSpecies];

    // Independent P(N_signal | activity, PDG).
    TH2D* signalMultiplicity = Ditto::SignalInjector::makeSignalMultiplicityHistogram(*tune, pdg, {-0.5, 0.5, 1.5});

    // P(pT | activity, PDG).
    std::vector<double> ptEdges;
    for (int i = 1; i <= pt_hists[iSpecies]->GetNbinsX(); ++i) {
      ptEdges.push_back(pt_hists[iSpecies]->GetXaxis()->GetBinLowEdge(i));
    }
    TH2D* signalPt = Ditto::SignalInjector::makePtHistogram(*tune, pdg, ptEdges);
    // P(eta | activity, PDG).
    std::vector<double> etaEdges;
    for (int i = 1; i <= eta_hists[iSpecies]->GetNbinsX(); ++i) {
      etaEdges.push_back(eta_hists[iSpecies]->GetXaxis()->GetBinLowEdge(i));
    }
    TH2D* signalEta = Ditto::SignalInjector::makeEtaHistogram(*tune, pdg, etaEdges);

    for (int iActivity = 1; iActivity <= nActivity; ++iActivity) {
      constexpr double injectionProbability = 0.10;
      signalMultiplicity->SetBinContent(iActivity, 1, 1.0 - injectionProbability); // N = 0
      signalMultiplicity->SetBinContent(iActivity, 2, injectionProbability);       // N = 1

      const double activityFraction = nActivity > 1 ? static_cast<double>(iActivity - 1) / static_cast<double>(nActivity - 1) : 0.0;

      for (int iPt = 1; iPt <= signalPt->GetNbinsY(); ++iPt) {
        signalPt->SetBinContent(iActivity, iPt, pt_hists[iSpecies]->GetBinContent(iPt));
      }

      for (int iEta = 1; iEta <= signalEta->GetNbinsY(); ++iEta) {
        signalEta->SetBinContent(iActivity, iEta, eta_hists[iSpecies]->GetBinContent(iEta));
      }
    }

    signalMultiplicity->Write();
    signalPt->Write();
    signalEta->Write();
  }

  output->Close();

  Printf("Wrote signal card to %s", outputFile);
}
