///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   Ditto.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Public interface for the Ditto generator.
///

#pragma once

#include <TClonesArray.h>
#include <TParticle.h>

#include <Pythia8/Pythia.h>

#include <cstdint>
#include <iosfwd>
#include <random>
#include <string>
#include <string_view>
#include <vector>

class TFile;
class TTree;

namespace Ditto
{

// Runtime options only.
//
// All physics configuration belongs to the tune file: beam IDs, energy,
// generated acceptance, particle species and the learned event/kinematic model.
struct Config {
  /// Ditto tune / generator-card file produced by Ditto::Tuner.
  std::string mTuneFile;

  /// Optional signal injection file. If empty, no signal is injected.
  std::string mSignalFile;

  /// Ditto random seed
  std::uint64_t mSeed = 1;

  // Measure total generator timing and print a summary when the Generator is
  // destroyed. Pure generation, TTree filling and PYTHIA export are reported
  // separately.
  bool mEnableTimingMetrics = false;

  // Add fine-grained timers inside generation. This is intended for profiling:
  // it performs many steady_clock::now() calls per particle and therefore
  // perturbs the absolute timing. Use enableTimingMetrics=true and this=false
  // for the least intrusive throughput benchmark.
  bool mEnableDetailedTimingMetrics = false;
};

struct EventInfo {
  std::uint64_t mEventNumber = 0;

  // Teacher beam configuration carried by the tune / generator card.
  int mBeamIdA = 0;      ///< PDG code of beam A, copied from the tune.
  int mBeamIdB = 0;      ///< PDG code of beam B, copied from the tune.
  double mSqrtSNN = 0.0; ///< Center-of-mass energy per nucleon pair, copied from the tune.

  // Exact event-level conditioning variables sampled from the tune.
  int mConditioningNch = -1;       ///< Exact number of charged particles in the activity acceptance, sampled from the tune.
  int mConditioningNSelected = -1; ///< Exact number of selected particles in the activity acceptance, sampled from the tune.
  int mActivityClass = -1;         ///< Exact activity class sampled from the tune, corresponding to the conditioning Nch.

  int mInjectedParticles = 0;     ///< Number of particles added by the optional signal overlay. The injected particles do not feed back into the background conditioning variables.
  int mInjectedDecayProducts = 0; ///< Number of decay products from the injected signal.
};

// Accumulated wall-clock timings. Values are stored in seconds.
//
// The detailed particle components are only populated when
// Config::enableDetailedTimingMetrics is true.
struct TimingMetrics {
  std::uint64_t mGeneratedEvents = 0;
  std::uint64_t mGeneratedParticles = 0;
  std::uint64_t mTreeFills = 0;
  std::uint64_t mPythiaExports = 0;

  // Coarse timings.
  double mGeneration = 0.0; // excludes TTree::Fill and PYTHIA export
  double mTreeFill = 0.0;
  double mPythiaExport = 0.0;

  // Event-level generation components.
  double mClearParticles = 0.0;
  double mEventInfo = 0.0;
  double mMultiplicity = 0.0;
  double mExpandArray = 0.0;
  double mParticleLoop = 0.0;

  // Particle-generation components.
  double mSpeciesSampling = 0.0;
  double mConstructedAt = 0.0;
  double mPtSampling = 0.0;
  double mEtaSampling = 0.0;
  double mPhiSampling = 0.0;
  double mMomentumMath = 0.0;
  double mParticleSetters = 0.0;

  // PYTHIA export components.
  double mPythiaReset = 0.0;
  double mPythiaAppend = 0.0;
  double mPythiaSystemSum = 0.0;
};

class Tune;
class SignalInjector;
struct TuneAliasSampler;
struct TuneDiscreteAliasSampler;
struct TuneRuntimeData;

class Generator
{
 public:
  explicit Generator(const Config& config);

  // Destructor, writes the TTree output if registered.
  ~Generator();

  // Generate and store one Ditto event internally.
  const EventInfo& generate();

  // Alias convenient in an event loop.
  const EventInfo& next() { return generate(); }

  const EventInfo& currentEvent() const { return mEvent; }
  const Config& config() const { return mConfig; }

  // Persistent tune object used as the complete physics generator card.
  const Tune& generatorCard() const;

  // Average pure generation time, excluding TTree filling and PYTHIA export.
  // Returns 0 if timing metrics are disabled or no event has been generated.
  double averageGenerationTimeUs() const;

  const TimingMetrics& timingMetrics() const { return mTiming; }
  void resetTimingMetrics();
  void printTimingMetrics(std::ostream& os) const;

  // Export the already-generated event. No new random numbers are consumed.
  // If reset=true, the target PYTHIA event is reset before filling.
  void loadParticles(Pythia8::Event& event, bool reset = true) const;

  /// @brief Creates an output file with a TTree containing the generated events
  /// @param fileName The name of the output file
  void registerTTreeOutput(const std::string& fileName);

 private:
  Config mConfig;
  EventInfo mEvent;
  std::mt19937_64 mRng;
  std::uint64_t mEventCounter = 0;

  mutable TimingMetrics mTiming;

  // The tune is mandatory and is the complete physics configuration.
  Tune* mTune = nullptr;                     //!
  TuneRuntimeData* mTuneRuntime = nullptr;   //!
  SignalInjector* mSignalInjector = nullptr; //!
  void validateConfig() const;
  void loadTune(const std::string& fileName);

  bool timingEnabled() const
  {
    return mConfig.mEnableTimingMetrics || mConfig.mEnableDetailedTimingMetrics;
  }

  bool detailedTimingEnabled() const
  {
    return mConfig.mEnableDetailedTimingMetrics;
  }

  double uniform01();
  double uniform(double min, double max);

  EventInfo makeEventInfo() const;
  int sampleMultiplicity(EventInfo& info);

  void sampleComposition(const EventInfo& info,
                         std::vector<int>& centralCounts,
                         std::vector<int>& otherCounts);

  double sampleAlias(const TuneAliasSampler& sampler);
  int sampleDiscreteAlias(const TuneDiscreteAliasSampler& sampler);
  double samplePhi();

  void makeParticle(TParticle* particle,
                    int pdg,
                    double mass,
                    double mass2,
                    double ptScale,
                    const TuneAliasSampler& ptSampler,
                    const TuneAliasSampler& etaSampler);

  // Output file for TTree writing. If empty, no output is written.
  static constexpr std::string_view kTreeName = "T";
  static constexpr std::string_view kBranchName = "Particles";
  TFile* mOutputFile = nullptr;
  TTree* mTree = nullptr;
  TClonesArray* mParticles = nullptr;
};

} // namespace Ditto
