HYPRE_LIB_DIR ?= $(HOME)/hypre/lib
HYPRE_INC_DIR ?=
FFTW_LIB_DIRS := $(patsubst %,-L%,$(subst :, ,$(LD_LIBRARY_PATH)))
FFTW_PKG_LIBS := $(shell pkg-config --libs fftw3 2>/dev/null)
FFTW_PKG_LIBS_SP := $(shell pkg-config --libs fftw3f 2>/dev/null)

override LIBS += -L$(HYPRE_LIB_DIR) -lHYPRE

ifneq ($(strip $(HYPRE_INC_DIR)),)
override INCS += -I$(HYPRE_INC_DIR)
endif

ifneq ($(strip $(FFTW_LIB_DIRS)),)
override LIBS += $(FFTW_LIB_DIRS)
endif

ifneq ($(strip $(FFTW_PKG_LIBS)),)
override LIBS += $(FFTW_PKG_LIBS)
else
override LIBS += -lfftw3
endif

ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3_threads
endif

ifeq ($(strip $(SINGLE_PRECISION)),1)
ifneq ($(strip $(FFTW_PKG_LIBS_SP)),)
override LIBS += $(FFTW_PKG_LIBS_SP)
else
override LIBS += -lfftw3f
endif
ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3f_threads
endif
endif
