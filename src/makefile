CC = /usr/bin/gcc
# CC = /usr/local/cuda/bin/nvcc
NVCC = /usr/local/cuda/bin/nvcc
GCCFLAGS = -c -Wall
NVCCFLAG = -L/usr/local/cuda/lib64 
ROOTDIR = .
SRCDIR = $(ROOTDIR)
OBJDIR = $(ROOTDIR)/obj
VPATH  = $(ROOTDIR):obj



OBJ = $(OBJDIR)/assess.o $(OBJDIR)/frequency.o $(OBJDIR)/blockFrequency.o \
      $(OBJDIR)/cusum.o $(OBJDIR)/runs.o $(OBJDIR)/longestRunOfOnes.o \
      $(OBJDIR)/serial.o $(OBJDIR)/rank.o $(OBJDIR)/discreteFourierTransform.o \
       $(OBJDIR)/universal.o $(OBJDIR)/Poker.o $(OBJDIR)/BinaryDerivative.o\
      $(OBJDIR)/approximateEntropy.o $(OBJDIR)/Autocorrelation.o\
       $(OBJDIR)/linearComplexity.o \
      $(OBJDIR)/dfft.o $(OBJDIR)/cephes.o $(OBJDIR)/matrix.o $(OBJDIR)/rats_app.o \
	  $(OBJDIR)/RunsDistribution.o $(OBJDIR)/linearComplexitycuda.o \
	  $(OBJDIR)/nonOverlappingTemplateMatchings.o
#  -L/usr/local/cuda/lib64 -lcuda -lcudart -lcublas -lcurand -lcudnn
assess: $(OBJ)
	$(CC) -o $@ $(OBJ) -lm -lcudart -lcuda $(NVCCFLAG) -lstdc++

$(OBJDIR)/assess.o: $(SRCDIR)/rats_demo.cpp rats.h rats_app.h
	$(CC) -o $@ -c $(SRCDIR)/rats_demo.cpp 

$(OBJDIR)/frequency.o: $(SRCDIR)/MonobitFrequency.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/MonobitFrequency.c

$(OBJDIR)/blockFrequency.o: $(SRCDIR)/BlockFrequency.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/BlockFrequency.c

$(OBJDIR)/cusum.o: $(SRCDIR)/Cumulative.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/Cumulative.c

$(OBJDIR)/runs.o: $(SRCDIR)/Runs.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/Runs.c

$(OBJDIR)/longestRunOfOnes.o: $(SRCDIR)/LongestRun.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/LongestRun.c

$(OBJDIR)/rank.o: $(SRCDIR)/MatrixRank.c cephes.h rats.h matrix.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/MatrixRank.c

$(OBJDIR)/discreteFourierTransform.o: $(SRCDIR)/DiscreteFourierTransform.c \
        cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/DiscreteFourierTransform.c

$(OBJDIR)/nonOverlappingTemplateMatchings.o: \
        $(SRCDIR)/nonOverlappingTemplateMatchings.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/nonOverlappingTemplateMatchings.c

# $(OBJDIR)/overlappingTemplateMatchings.o: \
#         $(SRCDIR)/overlappingTemplateMatchings.c defs.h externs.h utilities.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/overlappingTemplateMatchings.c

$(OBJDIR)/universal.o: $(SRCDIR)/MaurerUniversal.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/MaurerUniversal.c

$(OBJDIR)/approximateEntropy.o: $(SRCDIR)/ApproximateEntropy.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/ApproximateEntropy.c

# $(OBJDIR)/randomExcursions.o: $(SRCDIR)/randomExcursions.c defs.h externs.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/randomExcursions.c

# $(OBJDIR)/randomExcursionsVariant.o: $(SRCDIR)/randomExcursionsVariant.c defs.h externs.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/randomExcursionsVariant.c

$(OBJDIR)/serial.o: $(SRCDIR)/Serial.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/Serial.c

$(OBJDIR)/linearComplexity.o: $(SRCDIR)/LinearComplexity.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/LinearComplexity.c

$(OBJDIR)/linearComplexitycuda.o: $(SRCDIR)/LinearComplexityCUDA.cu  rats.h
	$(NVCC) -o $@ -c $(SRCDIR)/LinearComplexityCUDA.cu

$(OBJDIR)/dfft.o: $(SRCDIR)/dfft.c
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/dfft.c

$(OBJDIR)/matrix.o: $(SRCDIR)/matrix.c matrix.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/matrix.c

# $(OBJDIR)/genutils.o: $(SRCDIR)/genutils.c config.h genutils.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/genutils.c

$(OBJDIR)/cephes.o: $(SRCDIR)/cephes.c cephes.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/cephes.c

$(OBJDIR)/Poker.o: $(SRCDIR)/Poker.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/Poker.c

$(OBJDIR)/BinaryDerivative.o: $(SRCDIR)/BinaryDerivative.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/BinaryDerivative.c

$(OBJDIR)/Autocorrelation.o: $(SRCDIR)/Autocorrelation.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/Autocorrelation.c

$(OBJDIR)/rats_app.o: $(SRCDIR)/rats_app.c rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/rats_app.c

$(OBJDIR)/RunsDistribution.o: $(SRCDIR)/RunsDistribution.c cephes.h rats.h
	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/RunsDistribution.c
# $(OBJDIR)/utilities.o: $(SRCDIR)/utilities.c defs.h externs.h utilities.h config.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/utilities.c

# $(OBJDIR)/generators.o: $(SRCDIR)/generators.c defs.h externs.h utilities.h \
#         config.h generators.h
# 	$(CC) -o $@ $(GCCFLAGS) $(SRCDIR)/generators.c

clean:
	rm -f assess $(OBJDIR)/*.o

re: clean assess
