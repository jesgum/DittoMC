///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   example.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Example Ditto generation macro.
///

#include "Ditto.h"
#include "DittoTuneDownload.h"

#include <TClonesArray.h>
#include <TParticle.h>

#include <Pythia8/Pythia.h>

#include <chrono>
#include <iomanip>
#include <iostream>

void example(const int nEvents,
             const char* tuneFile,
             bool saveTTree)
{
  Ditto::Config cfg;

  // The tune is the complete physics generator card: beam IDs, energy,
  // acceptance, species and all learned distributions come from this file.
  cfg.mTuneFile = tuneFile;

  // Runtime-only Ditto options.
  cfg.mSeed = 12345;
  cfg.mEnableTimingMetrics = true;
  // Fine-grained profiling. Disable this for the least intrusive absolute benchmark.
  cfg.mEnableDetailedTimingMetrics = false;

  Ditto::Generator generator(cfg);

  // Register the output
  if (saveTTree) {
    generator.registerTTreeOutput("Generation/Ditto.root");
  }

  Pythia8::Pythia pythia;

  const auto startTime = std::chrono::steady_clock::now();

  for (int i = 0; i < nEvents; i++) {
    if (i % 1000 == 0) {
      const auto now = std::chrono::steady_clock::now();
      const double elapsedSec = std::chrono::duration<double>(now - startTime).count();
      const double avgSecPerEvent = (i > 0) ? elapsedSec / i : 0.0;
      const double etaSec = (i > 0) ? avgSecPerEvent * (nEvents - i) : 0.0;

      std::cout << "Event " << i << " / " << nEvents
                << " [elapsed: " << std::fixed << std::setprecision(1) << elapsedSec << "s"
                << ", avg: " << std::setprecision(4) << avgSecPerEvent * 1e3 << " ms/ev"
                << ", ETA: " << std::setprecision(1) << etaSec << "s]"
                << "\n";
    }
    generator.generate();
    if (i == nEvents - 1) {
      generator.loadParticles(pythia.event);
    }
  }

  const auto endTime = std::chrono::steady_clock::now();
  const double totalSec = std::chrono::duration<double>(endTime - startTime).count();
  std::cout << "Done: " << nEvents << " events in " << std::fixed << std::setprecision(2)
            << totalSec << "s (avg " << std::setprecision(4)
            << (totalSec / nEvents) * 1e3 << " ms/event)\n";

  std::cout << "PYTHIA entries: " << pythia.event.size() << "\n";
  pythia.event.list();
}

void example(int cfg = 0)
{
  int nEvents = 0;
  std::string tuneFile;
  switch (cfg) {
    case 0:
      nEvents = 1E8;
      tuneFile = Ditto::Tunes::resolve("https://github.com/njacazio/DittoMC/releases/download/v1.0.0/Ditto_tune_pythia8_inel_136tev.root");
      break;
    case 1:
      nEvents = 1E6;
      tuneFile = "https://github.com/njacazio/DittoMC/releases/download/v1.0.0/Ditto_tune_pythia8_OO_536.root";
      break;
    case 2:
      nEvents = 1E4;
      tuneFile = "https://github.com/njacazio/DittoMC/releases/download/v1.0.0/Ditto_tune_pythia8_PbPb_536tev.root";
      break;
    default:
      throw std::invalid_argument("example: unrecognized cfg value: " + std::to_string(cfg));
  }

  example(nEvents, tuneFile.c_str(), true);
}
