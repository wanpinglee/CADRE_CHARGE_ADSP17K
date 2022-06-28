
MASTER_DIR=$(shell pwd)
BIN_DIR=$(MASTER_DIR)/bin
LIB=$(MASTER_DIR)/libs

AUTOCONF = autoconf
AUTOHEADER = autoheader

CFLAGS:=
ifeq ($(mode), debug)
	CFLAGS:=$(CFLAGS) -O0 -g -DDEBUG
else
	CFLAGS:=$(CFLAGS) -O3
endif

CXXFLAGS:=-std=c++11 $(CFLAGS)
export $(CFLAGS)
export $(CXXFLAGS)

BCFTOOLS_DIR=$(LIB)/bcftools
BCFTOOLS=$(BIN_DIR)/bcftools
HTS_DIR=$(LIB)/htslib
HTS_LIB=$(HTS_DIR)/libhts.a

all: $(COMPACT) $(BCFTOOLS)
.PHONY: all

clean:
	@rm -rf $(BIN_DIR)
	$(MAKE) clean -C $(HTS_DIR)
	$(MAKE) clean -C $(BCFTOOLS_DIR)
.PHONY: clean

$(BCFTOOLS): $(HTS_LIB)
	@echo "- Building in bcftools"
	@cd $(BCFTOOLS_DIR) && $(AUTOHEADER) && $(AUTOCONF) && ./configure
	@$(MAKE) --no-print-directory -C $(BCFTOOLS_DIR)
	@mkdir -p $(BIN_DIR)
	@cp $(BCFTOOLS_DIR)/bcftools $@


$(HTS_LIB):
	@echo "- Building in htslib"
	@cd $(HTS_DIR) && $(AUTOHEADER) && $(AUTOCONF) && ./configure --disable-lzma --disable-lcurl
	@$(MAKE) --no-print-directory -C $(HTS_DIR)
