## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS   = fpdf mcpdf

fpdf_SOURCES    = fpdf.cc start_prob.cc

mcpdf_SOURCES   = mcpdf.cc start_prob.cc

INCLUDES = -I../include -I../.

DATE = $(shell date)

install-data-hook: 
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION
