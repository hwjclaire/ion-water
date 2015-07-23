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

OBJECTS=dipolemlwf.o  output-cccmdipole.o  output.o  output-spread.o sulfuric.o sulfate.o sulfuric+mlwf.o output-quad.o\
        output-ice.o output-ice192.o output-ice144.o output-ice128.o output-water.o

 LDFLAGS = ../lib/libbasic.a  $(LIBS)

 dipolemlwf.x: dipolemlwf.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output.x: output.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-quad.x: output-quad.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-cccmdipole.x: output-cccmdipole.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-dipole.x: output-dipole.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-vel.x: output-vel.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-pos.x: output-pos.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-O-WF.x: output-O-WF.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-water1.x: output-water1.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-spread.x: output-spread.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice.x: output-ice.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice-all.x: output-ice-all.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice128.x: output-ice128.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice144.x: output-ice144.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice144-list.x: output-ice144-list.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice144-quad-list.x: output-ice144-quad-list.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice144-polariz.x: output-ice144-polariz.o
	$(LD) $(DFLAGS) -o $@ $^ ../ewald/libewald.a $(LDFLAGS) 
 output-ice192.x: output-ice192.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice192-1.x: output-ice192-1.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice96.x: output-ice96.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-ice-bulk96.x: output-ice-bulk96.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-water.x: output-water.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-water-mlwf.x: output-water-mlwf.o
	$(LD) $(DFLAGS) -o $@ $^  ../ewald/libewald.a $(LDFLAGS)
 output-water-dist.x: output-water-dist.o
	$(LD) $(DFLAGS) -o $@ $^  ../ewald/libewald.a $(LDFLAGS)
 output-alumina.x: output-alumina.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-al2o3.x: output-al2o3.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-spread2.x: output-spread2.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 output-spread-forCui.x: output-spread-forCui.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 cui2andy.x: cui2andy.o
	$(LD) $(DFLAGS) -o $@ $^
 cui2coord.x: cui2coord.o
	$(LD) $(DFLAGS) -o $@ $^
 gro2andy.x: gro2andy.o
	$(LD) $(DFLAGS) -o $@ $^
 
 clean :
	rm -f *.o 
	rm -f *.x

