#===============================================================================
#   File: Makefile
#
#   Example Makefile for Pythia8 written for HF in STAR
#   Author: Thomas Ullrich
#   Last modified: September 9, 2008
#
#   This needs 3 packages installed: Pythia8, LHAPDF, and ROOT
#
#   This is setup for Mac OS X 10.5.4 but should be easy to 
#   adapt for other (Unix) platforms. In principle changing the
#   program name (PROGRAM), the location of PYTHIA (PYTHIAPATH),
#   and the LHAPDF libraries (LHAPDFPATH) should do the job.
#   Note that the environment variable ROOTSYS needs to be set.
#   Otherwise define it here in the makefile.
#===============================================================================
PROGRAM  =  NPEHDelPhiCorr
PYTHIAPATH   = /star/u/zbtang/myTools/pythia8142
#LHAPDFPATH   = /star/u/huangbc/package/local/pythia8/LHAPDF-6.1.4/lib
LHAPDFPATH   = /star/u/zbtang/myTools/lhapdf570
#ROOTSYS  = /star/u/zbtang/myTools/root522/
ROOTSYS  = /star/u/zbtang/myTools/root

CXX      =  g++
CXXFLAGS =  -m64 -fno-inline -O  -W -Wall
CPPFLAGS = -I$(PYTHIAPATH)/include -I$(ROOTSYS)/include
LDFLAGS  = -L$(PYTHIAPATH)/lib/archive -L$(ROOTSYS)/lib -L$(LHAPDFPATH)/lib -lLHAPDF -lpythia8 -llhapdfdummy -L$(ROOTSYS)/lib -lCore -lCint  -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lfreetype -lpthread -lm -ldl -lHist

$(PROGRAM):	$(PROGRAM).cpp Makefile
		$(CXX) $(CXXFLAGS) $(PROGRAM).cpp $(CPPFLAGS) $(LDFLAGS) -o $(PROGRAM) 

