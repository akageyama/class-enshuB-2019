
c_source_files := $(shell ls -1 *.c *.h)

source-dirs := $(dir $(c_source_files))
source-notdirs := $(notdir $(c_source_files))
source-suffix := $(suffix $(c_source_files))
source-basename := $(basename $(c_source_files))

test:
	@echo source-dirs = $(source-dirs)
	@echo source-notdirs = $(source-notdirs)
	@echo source-suffix = $(source-suffix)
	@echo source-basename = $(source-basename)
