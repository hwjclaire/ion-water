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
CC=CC
LD=$(CC)

EXEC=test.x
OBJECTS=Matrix.o test.o Context.o

CXXFLAGS= -O0 -DADD_
LDFLAGS=$(CXXFLAGS)
 $(EXEC):    $(OBJECTS)
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
 depend :
	makedepend -Y  *.[cCh]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f *.o 
	rm -f $(EXEC) 
