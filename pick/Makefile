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
CXX=g++
LD=$(CXX)

OBJECTS=pickHB.o

 LDFLAGS = ../lib/libbasic.a  $(LIBPATH) $(LIBS)

 pickHB.x: pickHB.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 pickWater.x: pickWater.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 pickHB-1.x: pickHB-1.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output.x: output.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-cccmdipole.x: output-cccmdipole.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-spread.x: output-spread.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-spread-forCui.x: output-spread-forCui.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 cui2andy.x: cui2andy.o
	$(LD) $(DFLAGS) -o $@ $^

 clean :
	rm -f *.o 
	rm -f *.x

