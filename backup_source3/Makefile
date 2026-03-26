CXX = g++
CXXFLAGS = -O3 -march=native -mavx2 -mfma -std=c++20
TARGET = elec_vdw_fixed.exe
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	del /Q $(OBJS) $(TARGET)
