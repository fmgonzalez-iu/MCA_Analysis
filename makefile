CC = g++
CFLAGS = `root-config --cflags` -std=c++11 -I/usr/include/mysql
LDFLAGS = `root-config --libs` -lRMySQL
objects := $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp))
#objects := $(wildcard ./src/*.cpp)

all: CFLAGS += -O3
all: analyzerForeach

debug: CFLAGS += -g
debug: analyzerForeach

analyzerForeach: AnalyzerForeach.cpp $(objects)
	$(CC) $(CFLAGS) -o AnalyzerForeach AnalyzerForeach.cpp $(objects) $(LDFLAGS)

clean:
	rm src/*.o

analyzer: Analyzer.cpp $(objects)
	$(CC) $(CFLAGS) -o analyzer Analyzer.cpp $(objects) $(LDFLAGS)

TTNE-Plotter: TTNE-Plotter.cpp $(objects)
	$(CC) $(CFLAGS) -o TTNE-Plotter TTNE-Plotter.cpp $(objects) $(LDFLAGS)

CoincExtractor: CoincExtractor.cpp $(objects)
	$(CC) $(CFLAGS) -o CoincExtractor CoincExtractor.cpp $(objects) $(LDFLAGS)
	
SinglesExtractor: SinglesExtractor.cpp $(objects)
	$(CC) $(CFLAGS) -o SinglesExtractor SinglesExtractor.cpp $(objects) $(LDFLAGS)

PhotonPlotter: PhotonPlotter.cpp $(objects)
	$(CC) $(CFLAGS) -o PhotonPlotter PhotonPlotter.cpp $(objects) $(LDFLAGS)

LifetimeAnalyzer: LifetimeAnalyzer.cpp $(objects)
	$(CC) $(CFLAGS) -o LifetimeAnalyzer LifetimeAnalyzer.cpp $(objects) $(LDFLAGS)	
	
AutomatedAnalyzer: AutomatedAnalyzer.cpp $(objects)
	$(CC) $(CFLAGS) -o AutomatedAnalyzer AutomatedAnalyzer.cpp $(objects) $(LDFLAGS)

src/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

chanDebugger: Channel_Debugger.cpp $(objects)
	$(CC) $(CFLAGS) -o Channel_Debugger Channel_Debugger.cpp $(objects) $(LDFLAGS)
