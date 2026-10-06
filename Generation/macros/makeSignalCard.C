/// Example signal-card writer for Ditto.
///
/// The signal multiplicity of every PDG is sampled independently.
/// The pT distribution is a thermal-like mT exponential:
///
///   dN/dpT ~ pT * exp[-(mT - m) / T]
///
/// with mT = sqrt(m^2 + pT^2).
///
/// The eta distribution uses a broad super-Gaussian, giving an approximately
/// flat central region and smooth tails, more representative of particle
/// production at LHC energies than a narrow Gaussian.
///
/// Run from a Ditto environment, for example:
///   root -l -b -q 'examples/makeSignalCard.C("tune.root","signal.root")'

#include "DittoSignal.h"
#include "DittoTune.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::vector<double> makeUniformEdges(double min, double max, int nBins)
{
  if (nBins <= 0 || max <= min) {
    throw std::invalid_argument("Invalid histogram binning");
  }

  std::vector<double> edges(static_cast<std::size_t>(nBins) + 1);
  const double width = (max - min) / static_cast<double>(nBins);

  for (int i = 0; i <= nBins; ++i) {
    edges[static_cast<std::size_t>(i)] = min + i * width;
  }

  return edges;
}

/// Thermal-like transverse-momentum spectrum.
///
/// If the invariant spectrum is approximately exponential in mT,
///
///   (1 / pT) dN/dpT ~ exp[-(mT - m) / T],
///
/// then
///
///   dN/dpT ~ pT * exp[-(mT - m) / T].
///
double ptWeight(double pt, double mass, double temperature)
{
  const double mt = std::hypot(pt, mass);
  return pt * std::exp(-(mt - mass) / temperature);
}

/// Broad midrapidity distribution with smooth tails.
///
/// exponent = 2 -> Gaussian
/// exponent > 2 -> increasingly flat central plateau.
///
/// exponent = 4 is a useful generic approximation for an LHC-like
/// pseudorapidity distribution over the central region.
double etaWeight(double eta, double width, double exponent = 4.0)
{
  const double x = std::abs(eta) / width;
  return std::exp(-0.5 * std::pow(x, exponent));
}

} // namespace

void makeSignalCard(const char* tuneFile, const char* outputFile)
{
  auto tune = Ditto::Tune::load(tuneFile);
  if (!tune || tune->mActivityEdges.size() < 2) {
    throw std::runtime_error("Could not load valid Ditto activity edges");
  }

  const auto& activityEdges = tune->mActivityEdges;
  const int nActivity = static_cast<int>(activityEdges.size()) - 1;

  // Example injected species.
  //
  // signalMass is also the registry of signal species: its X-axis labels
  // contain the corresponding PDG codes.
  const std::vector<int> pdgs = {421, -421};

  const std::vector<double> masses = {1.86484, 1.86484}; // GeV/c^2

  if (pdgs.size() != masses.size()) {
    throw std::runtime_error("PDG and mass vectors have different sizes");
  }

  std::unique_ptr<TFile> output(TFile::Open(outputFile, "RECREATE"));
  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("Could not create ") + outputFile);
  }

  // --------------------------------------------------------------------------
  // Species registry and masses
  // --------------------------------------------------------------------------

  TH1D* signalMass = Ditto::SignalInjector::makeSignalMassHistogram(*tune, pdgs, masses);

  signalMass->Write();

  // --------------------------------------------------------------------------
  // Kinematic binning
  // --------------------------------------------------------------------------

  // Fine enough that sampling inside the histogram gives a smooth distribution.
  const auto ptEdges = makeUniformEdges(0.0, 20.0, 400);
  const auto etaEdges = makeUniformEdges(-5.0, 5.0, 200);

  // --------------------------------------------------------------------------
  // Signal distributions
  // --------------------------------------------------------------------------

  for (std::size_t iSpecies = 0; iSpecies < pdgs.size(); ++iSpecies) {
    const int pdg = pdgs[iSpecies];
    const double mass = masses[iSpecies];

    // Independent P(N_signal | activity, PDG).
    //
    // This example allows either zero or one injected particle.
    TH2D* signalMultiplicity = Ditto::SignalInjector::makeSignalMultiplicityHistogram(*tune, pdg, {-0.5, 0.5, 1.5});

    // P(pT | activity, PDG).
    TH2D* signalPt = Ditto::SignalInjector::makePtHistogram(*tune, pdg, ptEdges);
    // P(eta | activity, PDG).
    TH2D* signalEta = Ditto::SignalInjector::makeEtaHistogram(*tune, pdg, etaEdges);

    for (int iActivity = 1; iActivity <= nActivity; ++iActivity) {

      // ----------------------------------------------------------------------
      // Multiplicity
      // ----------------------------------------------------------------------

      // Example: independently inject this species in 10% of events.
      //
      // Replace this with an activity-dependent probability if desired.
      constexpr double injectionProbability = 0.10;

      signalMultiplicity->SetBinContent(iActivity, 1, 1.0 - injectionProbability); // N = 0
      signalMultiplicity->SetBinContent(iActivity, 2, injectionProbability);       // N = 1

      // ----------------------------------------------------------------------
      // pT
      // ----------------------------------------------------------------------

      // Effective inverse slope in GeV.
      //
      // A mild activity dependence is included as an example: harder spectra
      // at larger activity. The independent variable is exactly the same
      // activity class used by Ditto itself.
      const double activityFraction = nActivity > 1 ? static_cast<double>(iActivity - 1) / static_cast<double>(nActivity - 1) : 0.0;

      const double temperature = 0.55 + 0.15 * activityFraction; // GeV

      for (int iPt = 1; iPt <= signalPt->GetNbinsY(); ++iPt) {
        const double pt = signalPt->GetYaxis()->GetBinCenter(iPt);
        const double binWidth = signalPt->GetYaxis()->GetBinWidth(iPt);
        // Histogram contents represent integrated bin weights.
        signalPt->SetBinContent(iActivity, iPt, ptWeight(pt, mass, temperature) * binWidth);
      }

      // ----------------------------------------------------------------------
      // eta
      // ----------------------------------------------------------------------

      // Broad central distribution.
      //
      // etaWidth controls approximately where the fall-off becomes important.
      // A very small activity dependence is included only as an example.
      const double etaWidth = 2.6 + 0.15 * activityFraction;
      for (int iEta = 1; iEta <= signalEta->GetNbinsY(); ++iEta) {
        const double eta = signalEta->GetYaxis()->GetBinCenter(iEta);
        const double binWidth = signalEta->GetYaxis()->GetBinWidth(iEta);
        signalEta->SetBinContent(iActivity, iEta, etaWeight(eta, etaWidth) * binWidth);
      }
    }

    signalMultiplicity->Write();
    signalPt->Write();
    signalEta->Write();
  }

  output->Close();

  Printf("Wrote signal card to %s", outputFile);
}

void makeSignalCard(const int mode)
{
  switch (mode) {
    case 0:
      makeSignalCard("Tuning/tunes/Ditto_tune_pythia8_inel_136tev.root", "signal.root");
      break;
    case 1:
      makeSignalCard("Tuning/tunes/Ditto_tune_pythia8_PbPb_536tev.root", "signal_pbpb.root");
      break;
    default:
      throw std::invalid_argument("Invalid mode for makeSignalCard");
  }
}
