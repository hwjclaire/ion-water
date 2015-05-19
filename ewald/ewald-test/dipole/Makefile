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
#include .mk
#------------------------------------------------------------------------------
#
CXX=g++
LD=g++

EXEC=testEwald.x
OBJECTS=Ewald.o testEwald.o 

CXXFLAGS= -Ofast -DADD_
LDFLAGS=$(CXXFLAGS) /opt/xt-libsci/11.0.06/gnu/46/istanbul/lib/libsci_gnu.a
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 depend :
	makedepend -Y  *.[cCh]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o 
	rm -f $(EXEC) 
#------------------------------------------------------------------------------
#  Cleanup object files
 tar :
	tar -czf water-src.tgz *.[Ch] Makefile *.mk

#------------------------------------------------------------------------------
# DO NOT DELETE

Ewald.o: D3vector.h Ewald.h
Ewald.o: D3vector.h
testEwald.o: D3vector.h Ewald.h
