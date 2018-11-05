#
# Makefile
# rgregoir, 2018-11-02 12:27
#

LIBS = -lcurl -lm -lz -lstdc++
INCLUDES = -L./third_party/libBigWig -I./third_party/libBigWig

all:
	clang --std=c++11 $(INCLUDES) main.cc third_party/libBigWig/libBigWig.a $(LIBS) -o bigwig-merge

debug:
	clang --std=c++11 $(INCLUDES) -g main.cc third_party/libBigWig/libBigWig.a $(LIBS) -o bigwig-merge

# vim:ft=make
