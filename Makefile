all: release
# Platform specifics
UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
  CXX_PLATFORM = -static-libgcc -static-libstdc++ 
else
  CXX_PLATFORM = -static
endif

EXEC = laser
UTIL_EXEC = 

DIR_EXEC = ./executable
DIR_EXEC_DBG = ./executable/dbg

$(DIR_EXEC):
	mkdir -p $@
$(DIR_EXEC_DBG):
	mkdir -p $@

##################################################
# Third-party libs.
TABIX_INC = third/tabix
TABIX_LIB = third/tabix/libtabix.a

EIGEN_INC = third/eigen
EIGEN_LIB =  # Eigen are header files only

PCRE_INC = third/pcre/include
PCRE_LIB = third/pcre/lib/libpcreposix.a third/pcre/lib/libpcre.a

GSL_INC = third/gsl/include
GSL_LIB = third/gsl/lib/libgsl.a

$(TABIX_INC) $(TABIX_LIB):
	(cd third; make tabix)
$(EIGEN_INC) $(EIGEN_LIB):
	(cd third; make eigen)
$(PCRE_INC) $(PCRE_LIB):
	(cd third; make pcre)
$(GSL_INC) $(GSL_LIB):
	(cd third; make gsl)

THIRD_INC = $(TABIX_INC) $(EIGEN_INC) $(PCRE_INC) $(GSL_INC)
THIRD_LIB = $(TABIX_LIB) $(PCRE_LIB) $(GSL_LIB)
##################################################
# Our libs.
BASE_INC = ./base
BASE_LIB = ./base/lib-base.a
BASE_LIB_DBG = ./base/lib-dbg-base.a

VCF_INC = ./libVcf
VCF_LIB = ./libVcf/lib-vcf.a
VCF_LIB_DBG = ./libVcf/lib-dbg-vcf.a

REGRESSION_INC = ./regression
REGRESSION_LIB = ./regression/lib-regression.a
REGRESSION_LIB_DBG = ./regression/lib-dbg-regression.a

GONCALO_INC = ./libsrc
GONCALO_LIB = ./libsrc/lib-goncalo.a
GONCALO_LIB_DBG = ./libsrc/lib-dbg-goncalo.a

$(BASE_LIB):
	(cd base; make)
$(BASE_LIB_DBG):
	(cd base; make debug)
$(VCF_LIB):
	(cd libVcf; make)
$(VCF_LIB_DBG):
	(cd libVcf; make debug)
$(REGRESSION_LIB): $(EIGEN_INC)
	(cd regression; make)
$(REGRESSION_LIB_DBG):
	(cd regression; make debug)
$(GONCALO_LIB):
	(cd libsrc; make)
$(GONCALO_LIB_DBG):
	(cd libsrc; make debug)

##################################################

INCLUDE = . $(THIRD_INC) $(BASE_INC)
LIB = $(BASE_LIB) $(THIRD_LIB)
LIB_DBG = $(BASE_LIB_DBG) $(THIRD_LIB)
CXX_INCLUDE = $(addprefix -I, $(INCLUDE))
CXX_LIB = $(LIB) -lz -lbz2 -lm
CXX_LIB_DBG = $(LIB_DBG) -lz -lbz2 -lm


DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug lib lib-dbg

# to build lib, we will use reverse order
# so third party lib is build first.
reverse = $(if $(1),$(call reverse,$(wordlist 2,$(words $(1)),$(1)))) $(firstword $(1))
lib: $(call reverse,$(LIB))
lib-dbg: $(call reverse,$(LIB_DBG))

release: CXX_FLAGS = -fopenmp -O2 -DNDEBUG $(DEFAULT_CXXFLAGS) $(CXX_PLATFORM)
release: $(DIR_EXEC)/$(EXEC) util
$(DIR_EXEC)/$(EXEC): lib \
                     Main.o \
                     ReadData.o \
                     |$(DIR_EXEC)
	g++ -o $@ Main.o ReadData.o $(CXX_FLAGS) $(CXX_LIB)

debug: CXX_FLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS) 
debug: $(DIR_EXEC_DBG)/$(EXEC) util-dbg
$(DIR_EXEC_DBG)/$(EXEC): lib-dbg \
                         Main.o \
                         ReadData.o \
                         | $(DIR_EXEC_DBG)
	g++ -o $@ Main.o ReadData.o $(CXX_FLAGS) $(CXX_LIB_DBG) 


##################################################
GitVersion.h: .git/HEAD .git/index
	echo "const char *gitVersion = \"$(shell git rev-parse HEAD)\";" > $@

-include Main.d
Main.o: Main.cpp GitVersion.h
	g++ -MMD -c $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__

-include ReadData.d
ReadData.o: ReadData.cpp GitVersion.h
	g++ -MMD -c $(CXX_FLAGS) $< $(CXX_INCLUDE) -D__ZLIB_AVAILABLE__


##################################################
# build utils
util: $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC))
define BUILD_util
  TAR := $(DIR_EXEC)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O2 $(DEFAULT_CXXFLAGS) $(CXX_PLATFORM)
  $$(TAR): $$(SRC) $(LIB) | $(DIR_EXEC)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util, $(s))))

util-dbg: $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))
define BUILD_util_dbg
  TAR := $(DIR_EXEC_DBG)/$(notdir $(basename $(1)))
  SRC := $(1).cpp
  -include  $(1).d
  $$(TAR): CXX_FLAGS = -O0 -ggdb $(DEFAULT_CXXFLAGS)
  $$(TAR): $$(SRC) $(LIB_DBG) | $(DIR_EXEC_DBG)
	g++ -MMD -o $$@ $$< $$(CXX_FLAGS) $(CXX_INCLUDE) $(CXX_LIB_DBG)
endef
$(foreach s, $(UTIL_EXEC), $(eval $(call BUILD_util_dbg, $(s))))


clean: 
	rm -rf *.o *.d $(EXEC) \
        $(addprefix $(DIR_EXEC)/,$(UTIL_EXEC)) \
        $(addprefix $(DIR_EXEC_DBG)/,$(UTIL_EXEC))

deepclean: clean
	(cd base; make clean)
	(cd regression; make clean)
	(cd libVcf; make clean)

test: test1 test2

testPCA.o: testPCA.cpp
	g++ -c $< $(CXX_FLAGS) $(CXX_INCLUDE)
test1: testPCA.o ReadData.o
	g++ -o $@ $^ $(CXX_FLAGS) $(CXX_LIB_DBG)

testProcrustes.o: testProcrustes.cpp
	g++ -c $< $(CXX_FLAGS) $(CXX_INCLUDE) 
test2: testProcrustes.o ReadData.o
	g++ -o $@ $^ $(CXX_FLAGS) $(CXX_LIB_DBG)

test3: rvtest
	./rvtest --inVcf test.vcf.gz --outVcf test3.vcf --peopleIncludeID P2,NotValid,P3 --peopleExcludeID P3
test4: rvtest
	./rvtest --inVcf test.vcf.gz --make-bed test.plink

# archive 
DATE=$(shell date '+%m%d')
tar:
	tar zvchf rvtest.$(DATE).tgz *.h Main.cpp tabix*tar.bz2 

# arg: Argument.h Argument.cpp
# 	g++ -g -o Argument Argument.cpp
# RangeList: RangeList_test.cpp RangeList.h RangeList_test_input
# 	g++ -c $(CXXFLAGS) RangeList_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ -lz
# 	g++ -o $@ RangeList_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm

# IO: IO_test.cpp IO.h 
# 	g++ -c $(CXXFLAGS) IO_test.cpp -I../statgen/lib/include -I. -D__ZLIB_AVAILABLE__ 
# 	g++ -o $@ IO_test.o $(TABIX_LIB) $(STATGEN_LIB)  -lz -lm -lbz2
README.md:README.wiki
doc: README.md
	java -jar third/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 
