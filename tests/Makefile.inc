#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#      Makefile for the unit tests included in the AENET package
#-----------------------------------------------------------------------
#+ This file is part of the AENET package.
#+
#+ Copyright (C) 2012-2016 Nongnuch Artrith and Alexander Urban
#+
#+ This program is free software: you can redistribute it and/or modify
#+ it under the terms of the GNU General Public License as published by
#+ the Free Software Foundation, either version 3 of the License, or
#+ (at your option) any later version.
#+
#+ This program is distributed in the hope that it will be useful, but
#+ WITHOUT ANY WARRANTY; without even the implied warranty of
#+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#+ General Public License for more details.
#+
#+ You should have received a copy of the GNU General Public License
#+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
# 2015-05-26 Alexander Urban (AU) and Nongnuch Artrith (NA)
#-----------------------------------------------------------------------

TST      = ./tests/
TSTBIN   = ../unittests/

TEST1    = $(TSTBIN)test_xsflib.x
TEST2    = $(TSTBIN)test_lclist.x
TEST3    = $(TSTBIN)test_geometry.x
TEST4    = $(TSTBIN)test_feedforward.x
TEST5    = $(TSTBIN)test_symmfunc.x
TEST6    = $(TSTBIN)test_sfsetup.x
TESTS    = $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5) $(TEST6)

#------------------------------- rules --------------------------------#

.PHONY: test tests

test  : tests
tests : $(TESTS)
	@for t in $(TESTS); do if [ -x $$t ]; then $$t; fi; done

%.o : $(TST)%.f90 unittest.o
	$(FC) $(FCFLAGS) $< -o $*.o

$(TSTBIN)%.x : %.o $(OBJECTS) unittest.o Makefile
	$(LD) $(LDFLAGS) -o $@ $< $(OBJECTS) unittest.o $(LIB) $(NUMLIB)
