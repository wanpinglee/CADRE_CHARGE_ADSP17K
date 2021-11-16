
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

COMPACT=$(BIN_DIR)/compact_vcf
BCFTOOLS_DIR=$(LIB)/bcftools-1.12
BCFTOOLS=$(BIN_DIR)/bcftools
HTS_DIR=$(LIB)/htslib-1.12
HTS_LIB=$(HTS_DIR)/libhts.a

all: $(COMPACT) $(BCFTOOLS)
.PHONY: all

clean:
	@rm -rf $(BIN_DIR)
	$(MAKE) clean -C $(LIB)/compact_vcf
	rm -rf $(HTS_DIR)
	rm -rf $(BCFTOOLS_DIR)
.PHONY: clean

$(COMPACT):
	@mkdir -p $(BIN_DIR)
	@echo "- Building in compact_vcf"
	@$(MAKE) --no-print-directory --directory=$(LIB)/compact_vcf
	@cp $(LIB)/compact_vcf/bin/compact_vcf $@

$(BCFTOOLS): $(HTS_LIB)
	@echo "- Building in bcftools"
	@cd $(LIB) && tar -zxvf $(LIB)/bcftools-1.12.tar.gz
	@cd $(BCFTOOLS_DIR) && $(AUTOHEADER) && $(AUTOCONF) && ./configure
	@$(MAKE) --no-print-directory -C $(BCFTOOLS_DIR)
	@cp $(BCFTOOLS_DIR)/bcftools $@


$(HTS_LIB):
	@echo "- Building in htslib"
	@cd $(LIB) && tar -zxvf $(LIB)/htslib-1.12.tar.gz
	@cd $(HTS_DIR) && $(AUTOHEADER) && $(AUTOCONF) && ./configure --disable-lzma --disable-lcurl
	@$(MAKE) --no-print-directory -C $(HTS_DIR)
