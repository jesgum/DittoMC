BUILD_DIR ?= build
CMAKE ?= cmake
JOBS ?=
BUILD_TYPE ?= Release
CXX_FLAGS_RELEASE ?= -O3 -DNDEBUG
PRODUCTION_CXX_FLAGS := -O3 -DNDEBUG -march=native -mtune=native

CMAKE_CONFIGURE_ARGS := -S . -B $(BUILD_DIR) \
	-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
	-DCMAKE_CXX_FLAGS_RELEASE="$(CXX_FLAGS_RELEASE)" \
	-DBUILD_TESTING=ON
ifneq ($(strip $(PYTHIA_ROOT)),)
CMAKE_CONFIGURE_ARGS += -DPYTHIA_ROOT=$(PYTHIA_ROOT)
endif

.PHONY: all configure build install production clean rebuild tune gen val valTune doc doc-open test tuneAO2D tunepp tuneOO tunePbPb tuneAll

all: build

configure:
	$(CMAKE) $(CMAKE_CONFIGURE_ARGS)

build: configure
	$(CMAKE) --build $(BUILD_DIR) --parallel $(if $(JOBS),$(JOBS),)

install: build
	$(CMAKE) --install $(BUILD_DIR)

test: build
	ctest --test-dir build --output-on-failure

production:
	$(MAKE) CXX_FLAGS_RELEASE="$(PRODUCTION_CXX_FLAGS)" install

clean:
	$(CMAKE) --build $(BUILD_DIR) --target clean
	rm -r lib
	rm -r include
	rm -r build

rebuild: clean install

tune: install
	root -l -b -q Tuning/exampleTuner.C

tunepp: install
	root -l -b -q 'Tuning/exampleTuner.C(0, 1E7)'

tuneOO: install
	root -l -b -q 'Tuning/exampleTuner.C(1, 1E6)'

tunePbPb: install
	root -l -b -q 'Tuning/exampleTuner.C(2, 1E5)'

tuneAll: tunepp tuneOO tunePbPb
	echo "All tuning completed."

tuneAO2D: install
	root -l -b -q Tuning/tuneAO2D.C

gen: install
	root -l -b -q Generation/example.C

val: install
	root -l -b -q Validation/ValidateDitto.C

valTune: install
	root -l -b -q ValidationTune/ValidateDittoTune.C

doc:
	doxygen Doxyfile

doc-open: doc
	xdg-open docs/generated/html/index.html

signal: install
	root -l -b -q Generation/macros/makeSignalCardFromDQ.C