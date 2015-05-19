#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
#
# This file is part of Qbox
#
# Qbox is distributed under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 2 of
# the License, or (at your option) any later version.
# See the file COPYING in the root directory of this distribution
# or <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.70 2010-05-10 20:51:00 fgygi Exp $
#------------------------------------------------------------------------------
#
#include kraken-gnu.mk
include ../$(TARGET).mk
#------------------------------------------------------------------------------
#

OBJECTS=dipolemlwf.o  output-cccmdipole.o  output.o  output-spread.o sulfuric.o sulfate.o sulfuric+mlwf.o \
       sulfuric-quad.o   sulfuric-quad-1st.o


DFLAGS = -DADD_ 
LDFLAGS = $(OPT) ../lib/libbasic.a $(LIBS)

 sulfate.x: sulfate.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 sulfate-spread.x: sulfate-spread.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 sulfuric.x: sulfuric.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 sulfuric+mlwf.x: sulfuric+mlwf.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 sulfuric-quad.x: sulfuric-quad.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 sulfuric-quad-1st.x: sulfuric-quad-1st.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 clean :
	rm -f *.o 
	rm -f *.x

