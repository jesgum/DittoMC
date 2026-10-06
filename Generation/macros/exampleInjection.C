///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   exampleInjection.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/29
/// \brief  Test macro for the Ditto signal-injection scheme.
///

#include "Ditto.h"
#include "DittoTuneDownload.h"

#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct SpeciesQA {
  int pdg = 0;

  TH2D* inputMultiplicity = nullptr;
  TH2D* inputPt = nullptr;
  TH2D* inputEta = nullptr;

  TH2D* generatedMultiplicity = nullptr;
  TH2D* generatedPt = nullptr;
  TH2D* generatedEta = nullptr;

  std::uint64_t injectedParticles = 0;
  std::uint64_t eventsWithInjection = 0;
};

std::string signalObjectName(const char* prefix, int pdg)
{
  return std::string(prefix) + std::to_string(pdg);
}

int pdgFromLabel(const TAxis& axis, int bin)
{
  const char* label = axis.GetBinLabel(bin);
  if (!label || !label[0]) {
    throw std::runtime_error("exampleInjection: empty PDG label in signalMass");
  }

  std::size_t consumed = 0;
  const std::string value(label);
  const int pdg = std::stoi(value, &consumed);

  if (consumed != value.size()) {
    throw std::runtime_error("exampleInjection: invalid PDG label '" + value + "'");
  }
  return pdg;
}

TH2D* cloneHistogram(const TH2D& source,
                     const std::string& name,
                     TFile& output,
                     bool reset)
{
  auto* histogram =
    static_cast<TH2D*>(source.Clone(name.c_str()));

  if (!histogram) {
    throw std::runtime_error("exampleInjection: could not clone histogram " + std::string(source.GetName()));
  }

  histogram->SetDirectory(&output);

  if (reset) {
    histogram->Reset("ICES");
  }

  return histogram;
}

} // namespace

void exampleInjection(const int nEvents,
                      const char* tuneFile,
                      const char* signalFile,
                      const char* outputFile,
                      const bool listLastEvent)
{
  if (nEvents <= 0) {
    throw std::invalid_argument("exampleInjection: nEvents must be positive");
  }

  // --------------------------------------------------------------------------
  // Read the signal card.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> signalInput(TFile::Open(signalFile, "READ"));
  if (!signalInput || signalInput->IsZombie()) {
    throw std::runtime_error(std::string("exampleInjection: could not open signal file ") + signalFile);
  }

  TH1D* signalMass = nullptr;
  signalInput->GetObject("signalMass", signalMass);

  if (!signalMass) {
    throw std::runtime_error("exampleInjection: signal file does not contain signalMass");
  }

  const int nSpecies = signalMass->GetNbinsX();
  if (nSpecies <= 0) {
    throw std::runtime_error("exampleInjection: signalMass contains no species");
  }

  // --------------------------------------------------------------------------
  // Configure Ditto exactly as in Generation/example.C, adding the signal card.
  // --------------------------------------------------------------------------

  Ditto::Config cfg;
  cfg.mTuneFile = tuneFile;
  cfg.mSignalFile = signalFile;
  cfg.mSeed = 12345;
  cfg.mEnableTimingMetrics = true;
  cfg.mEnableDetailedTimingMetrics = false;

  Ditto::Generator generator(cfg);

  // We use the PYTHIA export to inspect the generated event. Signal particles
  // are appended after the minimum-bias particles, so the last
  // EventInfo::mInjectedParticles entries are exactly the injected particles.
  Pythia8::Pythia pythia;

  // --------------------------------------------------------------------------
  // Prepare QA output.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> output(TFile::Open(outputFile, "RECREATE"));
  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("exampleInjection: could not create output file ") + outputFile);
  }

  std::vector<SpeciesQA> qa;
  qa.reserve(static_cast<std::size_t>(nSpecies));

  int maximumTotalMultiplicity = 0;

  for (int iSpecies = 1; iSpecies <= nSpecies; ++iSpecies) {
    SpeciesQA species;
    species.pdg = pdgFromLabel(*signalMass->GetXaxis(), iSpecies);

    const std::string multiplicityName = signalObjectName("signalMultiplicity_", species.pdg);
    const std::string ptName = signalObjectName("signalPt_", species.pdg);
    const std::string etaName = signalObjectName("signalEta_", species.pdg);

    TH2D* multiplicity = nullptr;
    TH2D* pt = nullptr;
    TH2D* eta = nullptr;

    signalInput->GetObject(multiplicityName.c_str(), multiplicity);
    signalInput->GetObject(ptName.c_str(), pt);
    signalInput->GetObject(etaName.c_str(), eta);

    if (!multiplicity || !pt || !eta) {
      throw std::runtime_error("exampleInjection: missing signal histograms for PDG " + std::to_string(species.pdg));
    }

    // Keep copies of the input distributions in the QA file.
    species.inputMultiplicity = cloneHistogram(*multiplicity, signalObjectName("inputMultiplicity_", species.pdg), *output, false);
    species.inputPt = cloneHistogram(*pt, signalObjectName("inputPt_", species.pdg), *output, false);
    species.inputEta = cloneHistogram(*eta, signalObjectName("inputEta_", species.pdg), *output, false);

    // Generated histograms have exactly the same axes as the card objects.
    species.generatedMultiplicity = cloneHistogram(*multiplicity, signalObjectName("generatedMultiplicity_", species.pdg), *output, true);
    species.generatedPt = cloneHistogram(*pt, signalObjectName("generatedPt_", species.pdg), *output, true);
    species.generatedEta = cloneHistogram(*eta, signalObjectName("generatedEta_", species.pdg), *output, true);

    const int lastMultiplicityBin = multiplicity->GetNbinsY();
    maximumTotalMultiplicity += std::max(0, static_cast<int>(std::llround(multiplicity->GetYaxis()->GetBinCenter(lastMultiplicityBin))));

    qa.push_back(species);
  }

  maximumTotalMultiplicity = std::max(1, maximumTotalMultiplicity);

  TH1D hInjectedParticles("generatedInjectedParticles", ";N_{injected};events", maximumTotalMultiplicity + 1, -0.5, static_cast<double>(maximumTotalMultiplicity) + 0.5);

  // Use the activity axis from the first per-PDG signal histogram.
  const TAxis* activityAxis = qa.front().inputMultiplicity->GetXaxis();
  std::vector<double> activityEdges(static_cast<std::size_t>(activityAxis->GetNbins()) + 1);
  for (int i = 1; i <= activityAxis->GetNbins(); ++i) {
    activityEdges[static_cast<std::size_t>(i - 1)] = activityAxis->GetBinLowEdge(i);
  }
  activityEdges.back() = activityAxis->GetBinUpEdge(activityAxis->GetNbins());
  std::vector<double> totalMultiplicityEdges(static_cast<std::size_t>(maximumTotalMultiplicity) + 2);

  for (int i = 0; i <= maximumTotalMultiplicity + 1; ++i) {
    totalMultiplicityEdges[static_cast<std::size_t>(i)] = static_cast<double>(i) - 0.5;
  }

  TH2D hInjectedParticlesVsActivity("generatedInjectedParticlesVsActivity", ";activity;N_{injected}", activityAxis->GetNbins(), activityEdges.data(), maximumTotalMultiplicity + 1, totalMultiplicityEdges.data());

  // --------------------------------------------------------------------------
  // Generate and test.
  // --------------------------------------------------------------------------

  std::uint64_t totalInjected = 0;
  std::uint64_t eventsWithAnySignal = 0;

  const auto startTime = std::chrono::steady_clock::now();

  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent % 10000 == 0) {
      const auto now = std::chrono::steady_clock::now();
      const double elapsedSec = std::chrono::duration<double>(now - startTime).count();
      const double avgSecPerEvent = iEvent > 0 ? elapsedSec / static_cast<double>(iEvent) : 0.0;
      const double etaSec = iEvent > 0 ? avgSecPerEvent * static_cast<double>(nEvents - iEvent) : 0.0;

      std::cout << "Event " << iEvent << " / " << nEvents
                << " [elapsed: " << std::fixed << std::setprecision(1)
                << elapsedSec << " s"
                << ", avg: " << std::setprecision(4)
                << avgSecPerEvent * 1.e3 << " ms/ev"
                << ", ETA: " << std::setprecision(1)
                << etaSec << " s]\n";
    }

    const Ditto::EventInfo& info = generator.generate();

    if (info.mActivityClass < 0 || info.mActivityClass >= activityAxis->GetNbins()) {
      throw std::runtime_error("exampleInjection: invalid activity class returned by Ditto");
    }

    generator.loadParticles(pythia.event);

    const int nInjected = info.mInjectedParticles;

    if (nInjected < 0) {
      throw std::runtime_error("exampleInjection: negative injected multiplicity");
    }

    // PYTHIA entry 0 is the system pseudo-particle. The real Ditto particles
    // follow, and the injected particles are the last nInjected entries.
    if (nInjected > pythia.event.size() - 1) {
      throw std::runtime_error("exampleInjection: injected multiplicity exceeds event size");
    }

    const int firstInjected = pythia.event.size() - nInjected;
    const double activityValue =
      activityAxis->GetBinCenter(info.mActivityClass + 1);

    std::vector<int> speciesMultiplicities(
      qa.size(),
      0);

    for (int iParticle = firstInjected; iParticle < pythia.event.size(); ++iParticle) {
      const auto& particle = pythia.event[iParticle];
      const int pdg = particle.id();

      const auto speciesIt = std::find_if(qa.begin(), qa.end(), [pdg](const SpeciesQA& species) {
        return species.pdg == pdg;
      });

      if (speciesIt == qa.end()) {
        throw std::runtime_error("exampleInjection: injected unexpected PDG " + std::to_string(pdg));
      }

      const std::size_t speciesIndex = static_cast<std::size_t>(std::distance(qa.begin(), speciesIt));
      ++speciesMultiplicities[speciesIndex];
      ++speciesIt->injectedParticles;

      speciesIt->generatedPt->Fill(activityValue, particle.pT());
      speciesIt->generatedEta->Fill(activityValue, particle.eta());
    }

    int reconstructedInjected = 0;

    for (std::size_t iSpecies = 0; iSpecies < qa.size(); ++iSpecies) {
      auto& species = qa[iSpecies];
      const int multiplicity = speciesMultiplicities[iSpecies];
      reconstructedInjected += multiplicity;
      species.generatedMultiplicity->Fill(activityValue, multiplicity);
      if (multiplicity > 0) {
        ++species.eventsWithInjection;
      }
    }

    if (reconstructedInjected != nInjected) {
      throw std::runtime_error("exampleInjection: sum of per-species multiplicities (" +
                               std::to_string(reconstructedInjected) + ") differs from EventInfo::mInjectedParticles (" + std::to_string(nInjected) + ")");
    }

    hInjectedParticles.Fill(nInjected);
    hInjectedParticlesVsActivity.Fill(activityValue, nInjected);
    totalInjected += static_cast<std::uint64_t>(nInjected);

    if (nInjected > 0) {
      ++eventsWithAnySignal;
    }
  }

  const auto endTime = std::chrono::steady_clock::now();
  const double totalSec = std::chrono::duration<double>(endTime - startTime).count();

  // --------------------------------------------------------------------------
  // Summary.
  // --------------------------------------------------------------------------

  std::cout << "\nSignal-injection test completed\n"
            << "  Events:             " << nEvents << "\n"
            << "  Injected particles: " << totalInjected << "\n"
            << "  <N injected>:       "
            << static_cast<double>(totalInjected) /
                 static_cast<double>(nEvents)
            << "\n"
            << "  Events with signal: "
            << eventsWithAnySignal << " / " << nEvents
            << " ("
            << 100.0 *
                 static_cast<double>(eventsWithAnySignal) /
                 static_cast<double>(nEvents)
            << "%)\n";

  for (const auto& species : qa) {
    std::cout << "  PDG " << species.pdg
              << ": particles = " << species.injectedParticles
              << ", <N> = "
              << static_cast<double>(species.injectedParticles) /
                   static_cast<double>(nEvents)
              << ", events with >=1 = "
              << species.eventsWithInjection
              << " ("
              << 100.0 *
                   static_cast<double>(species.eventsWithInjection) /
                   static_cast<double>(nEvents)
              << "%)\n";
  }

  std::cout << "  Runtime:            "
            << std::fixed << std::setprecision(2)
            << totalSec << " s\n"
            << "  QA output:          "
            << outputFile << "\n";

  if (listLastEvent) {
    std::cout << "\nLast generated event:\n";
    pythia.event.list();
  }

  output->cd();

  // Histograms cloned above belong to output. Explicit Write() makes the
  // intended contents clear in an interpreted ROOT macro.
  signalMass->Clone("inputSignalMass")->Write();

  for (const auto& species : qa) {
    species.inputMultiplicity->Write();
    species.inputPt->Write();
    species.inputEta->Write();

    species.generatedMultiplicity->Write();
    species.generatedPt->Write();
    species.generatedEta->Write();
  }

  hInjectedParticles.Write();
  hInjectedParticlesVsActivity.Write();

  output->Close();
}

// Convenience overload using the standard pp 13.6 TeV Ditto tune.
void exampleInjection(const int nEvents = 100000,
                      const char* signalFile = "https://github.com/njacazio/DittoMC/releases/download/v1.0.1/signal.root")
{
  const std::string tuneFile = Ditto::Tunes::resolve("Tuning/tunes/Ditto_tune_pythia8_inel_136tev.root");
  exampleInjection(nEvents, tuneFile.c_str(), signalFile, "Generation/DittoSignalTest.root", true);
}
