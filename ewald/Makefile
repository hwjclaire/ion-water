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

OBJECTS=Context.o  Ewald.o  Matrix.o  Polariz.o


 LDFLAGS = $(OPT) ../lib/libbasic.a $(LIBPATH) $(LIBS)


 EXEC=libewald.a
 $(EXEC):    $(OBJECTS)
	ar cr $@ $^ 
 testPolariz.x: testPolariz.o Polariz.o Matrix.o Context.o Ewald.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 testPolarizDM.x: testPolarizDM.o PolarizDM.o Matrix.o Context.o Ewald.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 NGPolariz.x: NGPolariz.o Polariz.o Matrix.o Context.o Ewald.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 clean :
	rm -f *.o 
	rm -f *.x

