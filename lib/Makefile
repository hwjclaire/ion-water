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

EXEC=libbasic.a
OBJECTS=NG.o  Stat.o  Water.o Sulfate.o D3vector.o Tensor.o Halide.o Alkali.o

 $(EXEC):    $(OBJECTS)
	ar cr $@ $^ 
 clean :
	rm -rf *.o *.a
