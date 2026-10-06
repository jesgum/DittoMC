///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoSignal.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/29
/// \brief  Internal configurable signal injector for Ditto.
///

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class TClonesArray;
class TH2D;
class TH1D;

namespace Ditto
{
class Tune;

// Internal helper used by Generator to overlay a configurable signal after the
// minimum-bias event has been generated.
//
// The signal card is conditioned on exactly the same activity-class axis used
// by the main Ditto kinematic model. Each signal species is injected
// independently, with its own multiplicity, pT and eta distributions. The
// injector owns a separate RNG stream so enabling signal injection does not
// perturb the minimum-bias RNG sequence.
class SignalInjector
{
 public:
  struct InjectionResult {
    int mInjectedParticles = 0;
    int mDecayProducts = 0;

    int totalParticles() const { return mInjectedParticles + mDecayProducts; }
  };

  SignalInjector(const std::string& fileName,
                 const std::vector<double>& activityEdges,
                 int finalStatus,
                 std::uint64_t seed);
  ~SignalInjector();

  SignalInjector(const SignalInjector&) = delete;
  SignalInjector& operator=(const SignalInjector&) = delete;
  SignalInjector(SignalInjector&&) = delete;
  SignalInjector& operator=(SignalInjector&&) = delete;

  // Append the injected particles starting from firstParticleIndex and return
  // the total number of particles added across all configured signal species.
  InjectionResult inject(TClonesArray& particles, int firstParticleIndex, int activityClass);

  // Helpers for constructing signal-card objects. The signalMass histogram is
  // the species registry: its X-axis labels are the PDG codes. Each registered
  // PDG must have signalMultiplicity_<PDG>, signalPt_<PDG> and signalEta_<PDG>.
  static TH2D* makeSignalMultiplicityHistogram(const Tune& tune, int pdg, const std::vector<double>& multiplicityEdges);
  static TH1D* makeSignalMassHistogram(const Tune& tune, const std::vector<int>& pdgs, const std::vector<double>& masses);
  static TH2D* makePtHistogram(const Tune& tune, int pdg, const std::vector<double>& ptEdges);
  static TH2D* makeEtaHistogram(const Tune& tune, int pdg, const std::vector<double>& etaEdges);

 private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace Ditto
