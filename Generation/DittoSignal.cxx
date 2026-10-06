///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoSignal.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/29
/// \brief  Implementation of the configurable Ditto signal injector.
///

#include "DittoSignal.h"

#include "DittoTune.h"

#include <TAxis.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TParticle.h>

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Ditto
{
namespace
{

constexpr double kTwoPi = 6.283185307179586476925286766559;
constexpr std::uint64_t kDecaySeedOffset = 0xd1b54a32d192ed03ULL;
constexpr std::uint64_t kMaxPythiaSeed = 900000000ULL;
constexpr const char* kMassName = "signalMass";
constexpr const char* kDecayName = "signalDecay";

std::string objectName(const char* prefix, int pdg)
{
  return std::string(prefix) + std::to_string(pdg);
}

int pythiaSeed(std::uint64_t seed)
{
  return static_cast<int>(seed % kMaxPythiaSeed) + 1;
}

bool nearlyEqual(double a, double b)
{
  const double scale = std::max({1.0, std::abs(a), std::abs(b)});
  return std::abs(a - b) <= 1.e-12 * scale;
}

bool axisMatchesEdges(const TAxis& axis,
                      const std::vector<double>& edges)
{
  if (edges.size() < 2 || axis.GetNbins() != static_cast<int>(edges.size()) - 1) {
    return false;
  }

  for (int i = 0; i < axis.GetNbins(); ++i) {
    if (!nearlyEqual(axis.GetBinLowEdge(i + 1), edges[static_cast<std::size_t>(i)])) {
      return false;
    }
  }
  return nearlyEqual(axis.GetBinUpEdge(axis.GetNbins()), edges.back());
}

int parsePdgLabel(const char* rawLabel,
                  const std::string& objectName,
                  int bin)
{
  const std::string label = rawLabel ? rawLabel : "";
  const auto first = label.find_first_not_of(" \t\n\r");
  const auto last = label.find_last_not_of(" \t\n\r");
  if (first == std::string::npos) {
    throw std::runtime_error("Ditto: empty PDG label in " + objectName + " bin " + std::to_string(bin));
  }

  const std::string trimmed = label.substr(first, last - first + 1);
  std::size_t consumed = 0;
  long long parsed = 0;
  try {
    parsed = std::stoll(trimmed, &consumed);
  } catch (const std::exception&) {
    throw std::runtime_error("Ditto: invalid PDG label '" + trimmed + "' in " + objectName);
  }

  if (consumed != trimmed.size() || parsed == 0 ||
      parsed < std::numeric_limits<int>::min() || parsed > std::numeric_limits<int>::max()) {
    throw std::runtime_error("Ditto: invalid PDG label '" + trimmed + "' in " + objectName);
  }
  return static_cast<int>(parsed);
}

struct ContinuousSampler {
  std::vector<double> lowEdges;
  std::vector<double> widths;
  std::discrete_distribution<std::size_t> binDistribution;

  bool empty() const { return lowEdges.empty(); }

  double sample(std::mt19937_64& rng)
  {
    if (empty()) {
      throw std::runtime_error("Ditto: attempted to sample an empty signal distribution");
    }
    const std::size_t bin = binDistribution(rng);
    const double u = static_cast<double>(rng() >> 11) * 0x1.0p-53;
    return lowEdges[bin] + widths[bin] * u;
  }
};

struct DiscreteSampler {
  std::vector<int> values;
  std::discrete_distribution<std::size_t> binDistribution;

  bool empty() const { return values.empty(); }

  int sample(std::mt19937_64& rng)
  {
    if (empty()) {
      throw std::runtime_error("Ditto: attempted to sample an empty signal distribution");
    }
    return values[binDistribution(rng)];
  }

  bool canSamplePositive() const
  {
    return std::any_of(values.begin(), values.end(), [](int value) {
      return value > 0;
    });
  }
};

ContinuousSampler makeContinuousSlice(
  const TH2D& histogram,
  int xBin,
  double minimum = -std::numeric_limits<double>::infinity(),
  double maximum = std::numeric_limits<double>::infinity())
{
  ContinuousSampler result;
  std::vector<double> weights;

  const TAxis* axis = histogram.GetYaxis();
  for (int i = 1; i <= axis->GetNbins(); ++i) {
    const double originalLow = axis->GetBinLowEdge(i);
    const double originalHigh = axis->GetBinUpEdge(i);
    const double originalWidth = originalHigh - originalLow;
    if (originalWidth <= 0.0) {
      continue;
    }

    const double low = std::max(originalLow, minimum);
    const double high = std::min(originalHigh, maximum);
    if (high <= low) {
      continue;
    }

    const double content = std::max(0.0, histogram.GetBinContent(xBin, i));
    if (content <= 0.0) {
      continue;
    }

    // Histogram contents are interpreted as integrated bin weights. If the
    // accepted range clips a boundary bin, keep the corresponding fraction.
    const double acceptedFraction = (high - low) / originalWidth;
    const double weight = content * acceptedFraction;
    if (weight <= 0.0) {
      continue;
    }

    result.lowEdges.push_back(low);
    result.widths.push_back(high - low);
    weights.push_back(weight);
  }

  if (!weights.empty()) {
    result.binDistribution = std::discrete_distribution<std::size_t>(weights.begin(), weights.end());
  }
  return result;
}

DiscreteSampler makeMultiplicitySlice(const TH2D& histogram,
                                      int xBin)
{
  DiscreteSampler result;
  std::vector<double> weights;

  const TAxis* axis = histogram.GetYaxis();
  for (int i = 1; i <= axis->GetNbins(); ++i) {
    const double weight = std::max(0.0, histogram.GetBinContent(xBin, i));
    if (weight <= 0.0) {
      continue;
    }

    const double value = axis->GetBinCenter(i);
    const long long rounded = std::llround(value);
    if (rounded < 0 || !nearlyEqual(value, static_cast<double>(rounded))) {
      throw std::runtime_error("Ditto: " + std::string(histogram.GetName()) +
                               " Y-bin centers with non-zero weight must be non-negative integers "
                               "(use e.g. edges -0.5, 0.5, 1.5, ...)");
    }
    if (rounded > std::numeric_limits<int>::max()) {
      throw std::runtime_error("Ditto: signal multiplicity is too large");
    }

    result.values.push_back(static_cast<int>(rounded));
    weights.push_back(weight);
  }

  if (!weights.empty()) {
    result.binDistribution = std::discrete_distribution<std::size_t>(weights.begin(), weights.end());
  }
  return result;
}

double uniformPhi(std::mt19937_64& rng)
{
  const double u = static_cast<double>(rng() >> 11) * 0x1.0p-53;
  return kTwoPi * u;
}

void fillParticle(TParticle& particle,
                  int pdg,
                  double mass,
                  double mass2,
                  int finalStatus,
                  double pt,
                  double eta,
                  double phi)
{
  double sinPhi = 0.0;
  double cosPhi = 0.0;
#if defined(__GNUC__) || defined(__clang__)
  __builtin_sincos(phi, &sinPhi, &cosPhi);
#else
  sinPhi = std::sin(phi);
  cosPhi = std::cos(phi);
#endif

  const double px = pt * cosPhi;
  const double py = pt * sinPhi;
  const double pz = pt * std::sinh(eta);
  const double p2 = pt * pt + pz * pz;
  const double energy = std::sqrt(p2 + mass2);

  particle.SetPdgCode(pdg);
  particle.SetStatusCode(finalStatus);
  particle.SetMother(0, -1);
  particle.SetMother(1, -1);
  particle.SetDaughter(0, -1);
  particle.SetDaughter(1, -1);
  particle.SetWeight(1.0f);
  particle.SetCalcMass(mass);
  particle.SetMomentum(px, py, pz, energy);
  particle.SetProductionVertex(0.0, 0.0, 0.0, 0.0);
  particle.SetPolarisation(0.0, 0.0, 0.0);
}

} // namespace

TH2D* SignalInjector::makeSignalMultiplicityHistogram(
  const Tune& tune,
  int pdg,
  const std::vector<double>& multiplicityEdges)
{
  if (multiplicityEdges.size() < 2) {
    throw std::invalid_argument("Ditto: signal multiplicity histogram requires at least two edges");
  }
  const auto& activityEdges = tune.mActivityEdges;
  if (activityEdges.size() < 2) {
    throw std::invalid_argument("Ditto: invalid tune activity edges");
  }
  const int nActivity = static_cast<int>(activityEdges.size()) - 1;
  const std::string name = objectName("signalMultiplicity_", pdg);
  const std::string title = ";N_{ch} activity;N_{signal} (PDG " + std::to_string(pdg) + ")";
  return new TH2D(name.c_str(), title.c_str(),
                  nActivity, activityEdges.data(),
                  static_cast<int>(multiplicityEdges.size()) - 1,
                  multiplicityEdges.data());
}

TH1D* SignalInjector::makeSignalMassHistogram(const Tune& tune,
                                              const std::vector<int>& pdgs,
                                              const std::vector<double>& masses)
{
  (void)tune;
  if (pdgs.empty() || pdgs.size() != masses.size()) {
    throw std::invalid_argument(
      "Ditto: signal PDG and mass lists must be non-empty and have the same size");
  }
  const int nSpecies = static_cast<int>(pdgs.size());
  TH1D* h = new TH1D(kMassName, ";signal species;mass (GeV/c^{2})",
                     nSpecies, 0.0, static_cast<double>(nSpecies));
  for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    const std::string label = std::to_string(pdgs[static_cast<std::size_t>(iSpecies)]);
    h->GetXaxis()->SetBinLabel(iSpecies + 1, label.c_str());
    h->SetBinContent(iSpecies + 1, masses[static_cast<std::size_t>(iSpecies)]);
  }
  return h;
}

TH2D* SignalInjector::makePtHistogram(const Tune& tune,
                                      int pdg,
                                      const std::vector<double>& ptEdges)
{
  if (ptEdges.size() < 2) {
    throw std::invalid_argument("Ditto: signal pT histogram requires at least two edges");
  }
  const auto& activityEdges = tune.mActivityEdges;
  if (activityEdges.size() < 2) {
    throw std::invalid_argument("Ditto: invalid tune activity edges");
  }
  const int nActivity = static_cast<int>(activityEdges.size()) - 1;
  const std::string name = objectName("signalPt_", pdg);
  const std::string title = ";N_{ch} activity;p_{T} (GeV/c) (PDG " + std::to_string(pdg) + ")";
  return new TH2D(name.c_str(), title.c_str(), nActivity, activityEdges.data(), static_cast<int>(ptEdges.size()) - 1, ptEdges.data());
}

TH2D* SignalInjector::makeEtaHistogram(const Tune& tune,
                                       int pdg,
                                       const std::vector<double>& etaEdges)
{
  if (etaEdges.size() < 2) {
    throw std::invalid_argument("Ditto: signal eta histogram requires at least two edges");
  }
  const auto& activityEdges = tune.mActivityEdges;
  if (activityEdges.size() < 2) {
    throw std::invalid_argument("Ditto: invalid tune activity edges");
  }
  const int nActivity = static_cast<int>(activityEdges.size()) - 1;
  const std::string name = objectName("signalEta_", pdg);
  const std::string title = ";N_{ch} activity;#eta (PDG " + std::to_string(pdg) + ")";
  return new TH2D(name.c_str(), title.c_str(), nActivity, activityEdges.data(), static_cast<int>(etaEdges.size()) - 1, etaEdges.data());
}

struct SignalInjector::Impl {
  struct SpeciesData {
    int pdg = 0;
    double mass = 0.0;
    double mass2 = 0.0;
    std::vector<DiscreteSampler> multiplicityGivenActivity;
    std::vector<ContinuousSampler> ptGivenActivity;
    std::vector<ContinuousSampler> etaGivenActivity;

    void print() const { Printf("Ditto: injected signal PDG %d, mass %.6f GeV/c^2\n", pdg, mass); }
  };

  explicit Impl(const std::string& fileName,
                const std::vector<double>& activityEdges,
                int status,
                std::uint64_t seed) : finalStatus(status), rng(seed)
  {
    if (activityEdges.size() < 2) {
      throw std::invalid_argument("Ditto: invalid activity edges for signal injection");
    }
    nActivityClasses = activityEdges.size() - 1;

    std::unique_ptr<TFile> input(TFile::Open(fileName.c_str(), "READ"));
    if (!input || input->IsZombie()) {
      throw std::runtime_error("Ditto: could not open signal file " + fileName);
    }

    TParameter<int>* decayParameter = nullptr;
    input->GetObject(kDecayName, decayParameter);
    if (decayParameter) {
      const int value = decayParameter->GetVal();
      if (value != 0 && value != 1) {
        throw std::runtime_error("Ditto: signalDecay must be either 0 or 1");
      }
      decaySignal = value != 0;
    }

    TH1D* massHistogram = nullptr;
    input->GetObject(kMassName, massHistogram);
    if (!massHistogram) {
      throw std::runtime_error("Ditto: signal file must contain TH1D 'signalMass'");
    }

    const int nSpecies = massHistogram->GetNbinsX();
    if (nSpecies <= 0) {
      throw std::runtime_error("Ditto: signalMass must contain at least one signal species");
    }

    std::unordered_set<int> seenPdgs;
    species.reserve(static_cast<std::size_t>(nSpecies));

    for (int iSpecies = 1; iSpecies <= nSpecies; ++iSpecies) {
      const int pdg = parsePdgLabel(massHistogram->GetXaxis()->GetBinLabel(iSpecies), kMassName, iSpecies);
      if (!seenPdgs.insert(pdg).second) {
        throw std::runtime_error("Ditto: duplicate signal PDG " + std::to_string(pdg));
      }

      const double mass = massHistogram->GetBinContent(iSpecies);
      if (!std::isfinite(mass) || mass < 0.0) {
        throw std::runtime_error("Ditto: invalid signal mass for PDG " + std::to_string(pdg));
      }

      const std::string multiplicityName = objectName("signalMultiplicity_", pdg);
      const std::string ptName = objectName("signalPt_", pdg);
      const std::string etaName = objectName("signalEta_", pdg);

      TH2D* multiplicityHistogram = nullptr;
      TH2D* ptHistogram = nullptr;
      TH2D* etaHistogram = nullptr;
      input->GetObject(multiplicityName.c_str(), multiplicityHistogram);
      input->GetObject(ptName.c_str(), ptHistogram);
      input->GetObject(etaName.c_str(), etaHistogram);

      if (!multiplicityHistogram || !ptHistogram || !etaHistogram) {
        throw std::runtime_error("Ditto: missing TH2D '" + multiplicityName + "', '" + ptName + "' or '" + etaName + "'");
      }

      if (!axisMatchesEdges(*multiplicityHistogram->GetXaxis(), activityEdges) ||
          !axisMatchesEdges(*ptHistogram->GetXaxis(), activityEdges) ||
          !axisMatchesEdges(*etaHistogram->GetXaxis(), activityEdges)) {
        throw std::runtime_error("Ditto: signal multiplicity/pT/eta X axes must exactly match the tune activity edges for PDG " +
                                 std::to_string(pdg));
      }

      SpeciesData data;
      data.pdg = pdg;
      data.mass = mass;
      data.mass2 = mass * mass;
      data.multiplicityGivenActivity.reserve(nActivityClasses);
      data.ptGivenActivity.reserve(nActivityClasses);
      data.etaGivenActivity.reserve(nActivityClasses);
      data.print();

      for (std::size_t iActivity = 0; iActivity < nActivityClasses; ++iActivity) {
        const int rootBin = static_cast<int>(iActivity) + 1;
        auto multiplicitySampler = makeMultiplicitySlice(*multiplicityHistogram, rootBin);
        if (multiplicitySampler.empty()) {
          throw std::runtime_error("Ditto: empty " + multiplicityName + " slice for activity class " + std::to_string(iActivity));
        }

        auto ptSampler = makeContinuousSlice(*ptHistogram, rootBin, 0.0);
        auto etaSampler = makeContinuousSlice(*etaHistogram, rootBin);
        if (multiplicitySampler.canSamplePositive() && (ptSampler.empty() || etaSampler.empty())) {
          throw std::runtime_error("Ditto: empty signal pT/eta distribution for PDG " +
                                   std::to_string(pdg) +
                                   " in an activity class where its "
                                   "multiplicity can be positive: " +
                                   std::to_string(iActivity));
        }

        data.multiplicityGivenActivity.push_back(std::move(multiplicitySampler));
        data.ptGivenActivity.push_back(std::move(ptSampler));
        data.etaGivenActivity.push_back(std::move(etaSampler));
      }

      species.push_back(std::move(data));
    }

    sampledMultiplicities.resize(species.size(), 0);

    if (decaySignal) {
      configureDecayer(seed + kDecaySeedOffset);
    }
  }

  void configureDecayer(std::uint64_t seed)
  {
    decayer = std::make_unique<Pythia8::Pythia>();

    if (!decayer->readString("ProcessLevel:all = off") ||
        !decayer->readString("Random:setSeed = on") ||
        !decayer->readString("Random:seed = " + std::to_string(pythiaSeed(seed))) ||
        !decayer->readString("Print:quiet = on")) {
      throw std::runtime_error("Ditto: could not configure PYTHIA signal decayer");
    }

    for (const auto& selected : species) {
      if (!decayer->particleData.isParticle(selected.pdg)) {
        throw std::runtime_error("Ditto: signal decay requested for unknown PYTHIA PDG " +
                                 std::to_string(selected.pdg));
      }
    }

    if (!decayer->init()) {
      throw std::runtime_error("Ditto: PYTHIA signal decayer initialization failed");
    }

    Printf("Ditto: signal decays enabled with PYTHIA\n");
  }

  int decay(TClonesArray& particles,
            int firstParticleIndex,
            int nParents)
  {
    if (!decayer || nParents <= 0) {
      return 0;
    }

    auto& event = decayer->event;
    event.reset();
    const int baseIndex = event.size();

    for (int i = 0; i < nParents; ++i) {
      const auto* particle = static_cast<const TParticle*>(particles.UncheckedAt(firstParticleIndex + i));
      if (!particle) {
        throw std::runtime_error("Ditto: missing injected particle before signal decay");
      }

      const int index = event.append(particle->GetPdgCode(),
                                     1,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     particle->Px(),
                                     particle->Py(),
                                     particle->Pz(),
                                     particle->Energy(),
                                     particle->GetCalcMass());
      event[index].vProd(particle->Vx(), particle->Vy(), particle->Vz(), particle->T());
    }

    if (!decayer->moreDecays()) {
      throw std::runtime_error("Ditto: PYTHIA signal decay failed");
    }

    const int nParticles = event.size() - baseIndex;
    if (nParticles < nParents) {
      throw std::runtime_error("Ditto: invalid PYTHIA signal decay record");
    }
    if (nParticles > std::numeric_limits<int>::max() - firstParticleIndex) {
      throw std::runtime_error("Ditto: signal decay particle index overflow");
    }

    const int requiredSize = firstParticleIndex + nParticles;
    if (requiredSize > particles.GetSize()) {
      particles.Expand(requiredSize);
    }

    const auto rootIndex = [baseIndex, firstParticleIndex](int pythiaIndex) {
      return pythiaIndex < baseIndex ? -1 : firstParticleIndex + pythiaIndex - baseIndex;
    };

    for (int i = baseIndex; i < event.size(); ++i) {
      const auto& source = event[i];
      auto* particle = static_cast<TParticle*>(particles.ConstructedAt(rootIndex(i)));

      particle->SetPdgCode(source.id());
      particle->SetStatusCode(source.isFinal() ? finalStatus : source.status());
      particle->SetMother(0, rootIndex(source.mother1()));
      particle->SetMother(1, rootIndex(source.mother2()));
      particle->SetDaughter(0, rootIndex(source.daughter1()));
      particle->SetDaughter(1, rootIndex(source.daughter2()));
      particle->SetWeight(1.0f);
      particle->SetCalcMass(source.m());
      particle->SetMomentum(source.px(), source.py(), source.pz(), source.e());
      particle->SetProductionVertex(source.xProd(), source.yProd(), source.zProd(), source.tProd());
      particle->SetPolarisation(0.0, 0.0, 0.0);
    }

    return nParticles - nParents;
  }

  SignalInjector::InjectionResult inject(TClonesArray& particles,
                                         int firstParticleIndex,
                                         int activityClass)
  {
    SignalInjector::InjectionResult result;
    if (firstParticleIndex < 0) {
      throw std::runtime_error("Ditto: invalid first particle index for signal injection");
    }
    if (activityClass < 0 || static_cast<std::size_t>(activityClass) >= nActivityClasses) {
      throw std::runtime_error("Ditto: invalid activity class for signal injection");
    }

    const std::size_t activity = static_cast<std::size_t>(activityClass);
    int totalSignal = 0;

    // Each configured species is an independent injection process, conditioned
    // only on the activity class selected by the main Ditto event.
    for (std::size_t iSpecies = 0; iSpecies < species.size(); ++iSpecies) {
      const int multiplicity = species[iSpecies].multiplicityGivenActivity[activity].sample(rng);
      sampledMultiplicities[iSpecies] = multiplicity;
      if (multiplicity > std::numeric_limits<int>::max() - totalSignal) {
        throw std::runtime_error("Ditto: total injected signal multiplicity is too large");
      }
      totalSignal += multiplicity;
    }

    if (totalSignal == 0) {
      return result;
    }
    if (totalSignal > std::numeric_limits<int>::max() - firstParticleIndex) {
      throw std::runtime_error("Ditto: injected particle index overflow");
    }

    const int requiredSize = firstParticleIndex + totalSignal;
    if (requiredSize > particles.GetSize()) {
      particles.Expand(requiredSize);
    }

    int particleIndex = firstParticleIndex;
    for (std::size_t iSpecies = 0; iSpecies < species.size(); ++iSpecies) {
      auto& selected = species[iSpecies];
      const int multiplicity = sampledMultiplicities[iSpecies];
      if (multiplicity == 0) {
        continue;
      }

      auto& ptSampler = selected.ptGivenActivity[activity];
      auto& etaSampler = selected.etaGivenActivity[activity];
      if (ptSampler.empty() || etaSampler.empty()) {
        throw std::runtime_error("Ditto: empty signal kinematic sampler for PDG " + std::to_string(selected.pdg));
      }

      for (int i = 0; i < multiplicity; ++i) {
        const double pt = ptSampler.sample(rng);
        const double eta = etaSampler.sample(rng);
        const double phi = uniformPhi(rng);

        auto* particle = static_cast<TParticle*>(particles.ConstructedAt(particleIndex++));
        fillParticle(*particle,
                     selected.pdg,
                     selected.mass,
                     selected.mass2,
                     finalStatus,
                     pt,
                     eta,
                     phi);
      }
    }

    result.mInjectedParticles = totalSignal;
    if (decaySignal) {
      result.mDecayProducts = decay(particles, firstParticleIndex, totalSignal);
    }
    return result;
  }

  int finalStatus = 1;
  std::mt19937_64 rng;
  bool decaySignal = false;
  std::unique_ptr<Pythia8::Pythia> decayer;
  std::size_t nActivityClasses = 0;
  std::vector<SpeciesData> species;
  std::vector<int> sampledMultiplicities;
};

SignalInjector::SignalInjector(const std::string& fileName,
                               const std::vector<double>& activityEdges,
                               int finalStatus,
                               std::uint64_t seed)
  : mImpl(std::make_unique<Impl>(fileName, activityEdges, finalStatus, seed))
{
  Printf("Ditto: signal injection initialized from file '%s'\n", fileName.c_str());
}

SignalInjector::~SignalInjector() = default;

SignalInjector::InjectionResult SignalInjector::inject(TClonesArray& particles,
                                                       int firstParticleIndex,
                                                       int activityClass)
{
  return mImpl->inject(particles, firstParticleIndex, activityClass);
}

} // namespace Ditto
